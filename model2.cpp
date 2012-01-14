
#include "model2.hpp"
#include "pointdev.hpp"
#include "cuba.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <cstring>
#include <cmath>
#include <cassert>
#include <alloca.h>
#include <fstream>
#include <limits>

using namespace std;

struct IntegralBounds;

typedef int (*IntegrandF) (const double xx[],double ff[], int ncomp, void *data);

ofstream iout("int.dat",ios::out);

struct IntegralBounds {
    double xmax;
    double xmin;
               
    double ymax;
    double ymin;
               
    double zmax;
    double zmin;

    void* data;
	IntegrandF integrand;
};




int cuhreIntegrand(const int *ndim, const double xx[],
					const int *ncomp, double ff[], void *voidbounds) {
    IntegralBounds* bounds = (IntegralBounds*) voidbounds;
    double xxprime[3];

	//Apply the transformation from the unit cube [0,1]^3. xp,yp,zp is
	//the dummy variable in molecular coordinates. The free variable
	//is intDet->xyz
	xxprime[0]=xx[0]*(bounds->xmax - bounds->xmin) + bounds->xmin;
	xxprime[1]=xx[1]*(bounds->ymax - bounds->ymin) + bounds->ymin;
	xxprime[2]=xx[2]*(bounds->zmax - bounds->zmin) + bounds->zmin;

	bounds->integrand(xxprime,ff,*ncomp,bounds->data);

    for(int i=0;i<*ncomp;i++){assert(isfinite(ff[i]));}

    return 0;
}


void cuhreIntegrate(IntegrandF f,IntegralBounds* bounds,unsigned long ncomp,double* integral,void* data) {
    const static int NDIM = 3;
    const static double EPSREL = 1e-2;
    const static double EPSABS = 0.00;
    const static int VERBOSE = 0;
    const static int LAST = 4;
    const static int MINEVAL = 1;
    const static int MAXEVAL = 5000000;
    const static int KEY = 11; //Use 11 point quadriture

	int nregions, neval, fail;
	double* error = (double*)alloca(ncomp * sizeof(double));
    double* prob  = (double*)alloca(ncomp * sizeof(double));

	bounds->integrand=f;
    bounds->data=data;

	Cuhre(NDIM, ncomp, cuhreIntegrand, (void*)bounds,
		  EPSREL, EPSABS, VERBOSE | LAST,
		  MINEVAL, MAXEVAL, KEY,
		  &nregions, &neval, &fail, integral, error, prob);

    double result_scaling = (bounds->xmax - bounds->xmin) *
		(bounds->ymax - bounds->ymin) *
		(bounds->zmax - bounds->zmin);
    for(unsigned long i = 0; i < ncomp; i++){integral[i]*=result_scaling;}

    for(unsigned long i = 0; i < ncomp;i++) {
        assert(isfinite(integral[i]));
    }
    if(fail != 0 || neval >= MAXEVAL)  {
        cerr << "Warning, an integral failed to converge, ncomp=" << ncomp << endl;
        for(unsigned long i = 0; i < ncomp; i++) {
            cout << "ff[" << i << "] = " << integral[i] << " +- " << (error[i]*result_scaling) 
                 << " " << (error[i]*result_scaling/integral[i])*100 << "%" << endl;
        }
    }
}


std::string name_param(int param) {
    switch(param) {
    case PARAM_X: return "x";
    case PARAM_Y: return "y";
    case PARAM_Z: return "z";
                 
    case PARAM_CHI1:  return "chi_1";
    case PARAM_CHI2:  return "chi_2";
    case PARAM_CHIXY: return "chi_xy";
    case PARAM_CHIXZ: return "chi_xz";
    case PARAM_CHIYZ: return "chi_yz";

    case PARAM_STDDEV: return "stddev";
    }
    assert(false);
    return "";
}

void numerical_derivative(Vec3d evalAt,const Model* model,const double* params,double* gradient) {
	double* params_mutable = (double*)alloca(model->size*sizeof(double));
	memcpy(params_mutable,params,model->size*sizeof(double));

    for(unsigned long i = 0;i<model->size;i++) {
		double h     = abs(params[i]*0.001);
		double result_plus,result_minus;


		params_mutable[i] = params[i] + h;
		model->modelf(evalAt,params_mutable,&result_plus ,NULL);

		params_mutable[i] = params[i] - h;
		model->modelf(evalAt,params_mutable,&result_minus,NULL);

		params_mutable[i] = params[i];

		gradient[i] = (result_plus-result_minus)/(2*h);
	}
}


void eval_point(Vec3d evalAt,const double* pm,double* value, double* gradient) {
    double x = evalAt.x() - pm[PARAM_X];
    double y = evalAt.y() - pm[PARAM_Y];
    double z = evalAt.z() - pm[PARAM_Z];

	if(x == 0 && y == 0 && z == 0) {
		*value = 0;
        if(gradient != NULL) {
            for(unsigned long i = 0;i < 8; i++){gradient[i] = 0;}
        }
		return;
	}

    double chi_1 =  pm[PARAM_CHI1]; 
    double chi_2 =  pm[PARAM_CHI2];
    double chi_xy = pm[PARAM_CHIXY];
    double chi_xz = pm[PARAM_CHIXZ];
    double chi_yz = pm[PARAM_CHIYZ];

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;

    double xy = x*y;
    double xz = x*z;
    double yz = y*z;

    double r2 = x2+y2+z2;

    double r = sqrt(r2);
    double r5 = r2*r2*r;

    double inv12PiR5 = 1/(12*M_PI*r5);

    double A = (r2-3*x2)*chi_1 + (z2-y2)*chi_2 + 6*(xy*chi_xy + xz*chi_xz + yz*chi_yz);
    *value = inv12PiR5 * A; assert(isfinite(*value));

    if(gradient != NULL) {
        gradient[PARAM_CHI1]  = inv12PiR5*(r2-3*x2);
        gradient[PARAM_CHI2]  = inv12PiR5*(z2-y2);
        gradient[PARAM_CHIXY] = 6*inv12PiR5*xy;
        gradient[PARAM_CHIXZ] = 6*inv12PiR5*xz;
        gradient[PARAM_CHIYZ] = 6*inv12PiR5*yz;


        double fiveAOver12Pi7 = 5*A/(12*M_PI*r5*r2);

        //There is a minus here because this is a deriative with
        //respect to the position of the metal, which appears as
        //sigma(r - r_m) where sigma is a model assuming the metal
        //sits at the origin and r_m is the location of the metal
        gradient[PARAM_X] = -(-x*fiveAOver12Pi7 + (-4*x*chi_1          + 6*(y*chi_xy + z*chi_xz))*inv12PiR5);
        gradient[PARAM_Y] = -(-y*fiveAOver12Pi7 + (2*y*(chi_1 - chi_2) + 6*(x*chi_xy + z*chi_yz))*inv12PiR5);
        gradient[PARAM_Z] = -(-z*fiveAOver12Pi7 + (2*z*(chi_1 + chi_2) + 6*(x*chi_xz + y*chi_yz))*inv12PiR5);
	
        for(unsigned long i = 0;i<8;i++) {
            assert(isfinite(gradient[i]));
            assert(isfinite(gradient[i]*gradient[i]));
        }
    }
}

struct Userdata {
    //For sending cuhre.
    const Model* point_model;
    const double* params; 
	Vec3d evalAt;

    double reg_radius2;
};




int parse_params_file(const std::string& filename,const Model** model,std::vector<double>* params) {
    ifstream fin(filename.c_str());
    if(!fin.is_open()) {
        return PARAM_FILE_NOT_FOUND;
    }
    string model_name;
    fin >> model_name;
    if(model_name == "point") {
        *model = &point_model;
    } else if(model_name == "gauss_test") {
        cout << "Using gauss_test model " << endl;
        *model = &gaussian_model_testing;
    } else if(model_name == "gauss") {
        *model = &gaussian_model;
    } else {
        return  UNKNOWN_MODEL;
    }
    params->resize((*model)->size);
    for(unsigned long i = 0;i < (*model)->size; i++) {
        fin >> params->at(i);
        if(fin.eof()) {
            return NOT_ENOUGH_PARAMS;
        }
    }
    return PARSE_SUCESS;
}


void random_data(PRNG& prng,const Model& model,const double* params,unsigned long natoms,Dataset* dataset) {
	RandomDist rand;

    //We should center our imaginary spins around the dipole or
    //fitting will be very hard.
    double x = params[PARAM_X];
    double y = params[PARAM_Y];
    double z = params[PARAM_Z];
    
    dataset->nuclei.resize(natoms);
    dataset->vals.resize(natoms);

    for(unsigned long i = 0; i < natoms;i++) {
		dataset->nuclei[i] = Vec3d(x + 10*rand(prng),y + 10*rand(prng),z + 10*rand(prng));
        model.modelf(dataset->nuclei[i],params,&(dataset->vals[i]),NULL);
    }
}

//================================================================================//

int Integrand(const double xx[],double ff[],int ncomp, void* void_userdata) {
    Userdata* userdata = (Userdata*)(void_userdata);

    const double* params = userdata->params;
    double stddev = abs(params[PARAM_STDDEV]);

    unsigned long point_size = userdata->point_model->size;

    Vec3d xyz = Vec3d(xx);
	double r2 = xyz.r2();

    Vec3d metal = Vec3d(params);

    //The location of the singularity (optimisation: take these
    //calculations out of the integral)
    Vec3d singularity = userdata->evalAt - metal;
    double singularity_r2 = singularity.r2();

    //Potential optimisation: pull the calculation of these out of the integral
    double a_coef = 1/(stddev*stddev);   assert(isfinite(a_coef));
    double normalizer = pow(M_PI*stddev*stddev,-1.5);
    double theExp0 = exp(-a_coef*(singularity_r2));
    double rho0 = normalizer*theExp0;

    double stddev2 = stddev*stddev;

    //First derivative
    Vec3d rho0_grad = singularity *(-2/stddev2 * rho0);

    //A vector pointing from the centre of 1/r^3 to the center of the gaussian
    Vec3d expand_around = xyz - singularity;
    double expand_r2 = expand_around.r2();

    //Evauate the model
    double f;
    double* gradient = ncomp == 1 ? NULL: (double*)alloca(point_size*sizeof(double));

    Vec3d x_minus_xprime = userdata->evalAt - xyz;

    userdata->point_model->modelf(x_minus_xprime,params,&f,gradient);


    double theExp = exp(-a_coef*r2);
    double rho = normalizer*theExp;

    //Is any sort of regualisation needed?
    bool close      = expand_r2 < 4*userdata->reg_radius2;
    bool very_close = expand_r2 < userdata->reg_radius2;

    double correction = 0;

    if(close) {
        if(very_close) {
            correction = 1;
        } else {
            double reg_radius = sqrt(userdata->reg_radius2);
            double expand_r = sqrt(expand_r2);

            double t =  (expand_r - reg_radius)/(reg_radius);

            correction = magic_polynomial(1-t);
            //cout << 1-t << " " << correction << endl;;
        }
    }

    assert(isfinite(rho));
    assert(isfinite(rho0));
    assert(isfinite(f));

    double toSub =  correction*(rho0 + rho0_grad.dot(expand_around));

    if( abs((rho-toSub)/rho) < 1e-10) {
        //cout << " Warning, possible loss of precision, rho = "
        //<< rho << " toSub = " << toSub << "  diff = " << (rho-toSub) << endl;

        //cout << " computed diff = " << normalizer*gaussian_error_term_one(r,a,stddev) << endl;
    }
    ff[0] = (rho-toSub)*f;

    assert(isfinite(ff[0]));
    if(ncomp == 1) {
        //If we don't need a gradient, we can stop here.
        return 0;
    }

    for(unsigned long i=1;i<9;i++){
        ff[i] = (rho-toSub) * gradient[i-1];

        //if(i == 1) cout << sqrt(prime_to_singularity2) << ", " << ff[i] << endl;
    }
    ff[1] = 0;
    ff[2] = 0;
    ff[3] = 0;
    for(int i=0;i<ncomp;i++){assert(isfinite(ff[i]));}
    return 0;
}


void eval_gaussian(Vec3d evalAt,const double* params,double* value, double* gradient) {
    Userdata userdata;
    userdata.point_model = &point_model;
    userdata.params = params;
	userdata.evalAt = evalAt;

    double stddev = params[PARAM_STDDEV];

    IntegralBounds bounds;
    bounds.xmax =  7*stddev;
    bounds.xmin = -7*stddev;
    
    bounds.ymax =  7*stddev;
    bounds.ymin = -7*stddev;
    
    bounds.zmax =  7*stddev;
    bounds.zmin = -7*stddev;

    //Make sure that the region where regualization is applied is
    //wholey included or wholey excluded.

    double reg_radius    = stddev/3.0;
    userdata.reg_radius2 = reg_radius*reg_radius;
    unsigned long ncomp = gradient == NULL ? 1 : (POINT_SIZE+1);

	double* integral = (double*)alloca(ncomp*sizeof(double));

	cuhreIntegrate(Integrand,&bounds,ncomp,integral,(void*)&userdata);

    *value = integral[0];

    if(gradient != NULL) {
        memcpy(gradient,integral+1,POINT_SIZE*sizeof(double));

        double h = 0.0002;
        double val_p,val_m;
        double val_2p,val_2m;
        double val_3p,val_3m;
        double val_4p,val_4m;
        
        double mute_params[GAUSS_SIZE];
        userdata.params = mute_params;

        memcpy(mute_params,params,sizeof(double)*GAUSS_SIZE);

        mute_params[PARAM_STDDEV] = stddev + 4*h;
        cuhreIntegrate(Integrand,&bounds,1,&val_4p,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev + 3*h;
        cuhreIntegrate(Integrand,&bounds,1,&val_3p,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev + 2*h;
        cuhreIntegrate(Integrand,&bounds,1,&val_2p,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev + h;
        cuhreIntegrate(Integrand,&bounds,1,&val_p,(void*)&userdata);


        mute_params[PARAM_STDDEV] = stddev - h;
        cuhreIntegrate(Integrand,&bounds,1,&val_m,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev - 2*h;
        cuhreIntegrate(Integrand,&bounds,1,&val_2m,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev - 3*h;
        cuhreIntegrate(Integrand,&bounds,1,&val_3m,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev - 4*h;
        cuhreIntegrate(Integrand,&bounds,1,&val_4m,(void*)&userdata);

        //double fd1gradient = (val_p-val_m)/(2*h);
        //double fd2gradient = (val_2m - 8*val_m + 8*val_p - val_2p)/(12*h);
        //double fd3gradient = (-val_3m + 9*val_2m - 45*val_m + 45*val_p - 9*val_2p + val_3p)/(60*h);
        double fd4gradient = ( val_4m/280 - 4*val_3m/105 + val_2m/5 - 4*val_m/5
                              -val_4p/280 + 4*val_3p/105 - val_2p/5 + 4*val_p/5)/h;

        /*cout << "sgradient = " << sgradient
             << " fd1gradient = " << fd1gradient
             << " fd2gradient = " << fd2gradient
             << " fd3gradient = " << fd3gradient
             << " fd4gradient = " << fd4gradient
             << endl;*/

        gradient[PARAM_STDDEV] = fd4gradient;
    }
}

//================================================================================//

void eval_gaussian_testing(Vec3d evalAt,const double* params,double* value, double* gradient) {
    if(gradient != NULL) {
        //Don't bother with an analytical gradient
        for(unsigned long i = 0; i < gaussian_model_testing.size;i++) {
            gradient[i] = 0;
        }
    }

    double normalizer = 0;
    double total      = 0;

    double stddev = params[PARAM_STDDEV];
    double a_coef = 1/(stddev*stddev);  
    /*
      Eval (2l+1)^3 times on a squard grid and spacings of stddev/4
     */

    long l = 5;

    for(long i = -l; i < l; i++) { // -5,-4, ... , 4,5 if l=5
        for(long j = -l; j < l; j++) {
            for(long k = -l; k < l; k++) {
                double f;

                double x = i * (stddev/2);
                double y = j * (stddev/2);
                double z = k * (stddev/2);
                
                double r = sqrt(x*x+y*y+z*z);
                
                if(r < stddev/4) {
                    continue;
                }

                double rho = exp(-a_coef*r*r);

                Vec3d x_minus_xprime(evalAt.x() - x,
                                     evalAt.y() - y,
                                     evalAt.z() - z);
                eval_point(x_minus_xprime,params,&f,NULL);

                total += f * rho;
                normalizer += rho;
            }
        }
    }
    *value = total / normalizer;
}


//================================================================================//
int IntegrandNumDev(const double xx[],double ff[],int ncomp, void* void_userdata) {
    Userdata* userdata = (Userdata*)(void_userdata);
    Userdata userdata_copy = *userdata;

    unsigned long size = (userdata->point_model->size+1);

	double* params_mutable = (double*)alloca(size*sizeof(double));
	memcpy(params_mutable,userdata->params,size*sizeof(double));
    userdata_copy.params=params_mutable;

    assert(ncomp == int(1 + userdata->point_model->size + 1));

    Integrand(xx,ff,1,(void*)&userdata_copy);

    for(unsigned long i = 0; i < size; i++) {
        double param = userdata->params[i];

		double h     = abs(param*0.00001);
		double result_plus  = numeric_limits<double>::quiet_NaN();
        double result_minus = numeric_limits<double>::quiet_NaN();

		params_mutable[i] = param + h;
        Integrand(xx,&result_plus,1,(void*)&userdata_copy);

		params_mutable[i] = param - h;
        Integrand(xx,&result_minus,1,(void*)&userdata_copy);

		params_mutable[i] = param;

		ff[i+1] = (result_plus-result_minus)/(2*h);
        //if(i == 8 && abs(xx[2]) < 0.01) cout << xx[0] << " " << xx[1] << " " <<  ff[i+1] << endl;
	}

    return 0;
}


void eval_gaussian_num_dev(Vec3d evalAt,const double* params,double* value, double* gradient) {
    Userdata userdata;
    userdata.point_model = &point_model;
    userdata.params = params;
	userdata.evalAt = evalAt;
    double stddev = params[PARAM_STDDEV];

    IntegralBounds bounds;
    bounds.xmax =   5*stddev;
    bounds.xmin =  -5*stddev;
    
    bounds.ymax =   5*stddev;
    bounds.ymin =  -5*stddev;
    
    bounds.zmax =   5*stddev;
    bounds.zmin =  -5*stddev;


    double reg_radius    = stddev/2;
    userdata.reg_radius2 = reg_radius*reg_radius;

    //We we're doing the gradient, how many functions do we need? 
    unsigned long ncomp = gradient == NULL ? 1 : 1 + userdata.point_model->size + 1;

	double* integral = (double*)alloca(ncomp*sizeof(double));

    if(gradient == NULL) {
        cuhreIntegrate(Integrand      ,&bounds,ncomp,integral,(void*)&userdata);
    } else {
        cuhreIntegrate(IntegrandNumDev,&bounds,ncomp,integral,(void*)&userdata);
    }

    *value = integral[0];

    if(gradient != NULL) {
        memcpy(gradient,integral+1,8*sizeof(double));

        double s5 = stddev*stddev*stddev*stddev*stddev;

        double sigma_xxxx = eval_point_model_dev_xyz(params,2,0,0,evalAt);
        double sigma_yyyy = eval_point_model_dev_xyz(params,0,2,0,evalAt);
        double sigma_zzzz = eval_point_model_dev_xyz(params,0,0,2,evalAt);

        double sigma_xxyy = eval_point_model_dev_xyz(params,2,0,0,evalAt);
        double sigma_xxzz = eval_point_model_dev_xyz(params,0,2,0,evalAt);
        double sigma_yyzz = eval_point_model_dev_xyz(params,0,0,2,evalAt);

        gradient[PARAM_STDDEV] = 6*s5*(
                                       3.0/4*(sigma_xxxx + sigma_yyyy + sigma_zzzz)+
                                       1.0/2*(sigma_xxyy + sigma_xxzz + sigma_yyzz)
                                       )/24;
    }
}

const Model point_model    = {eval_point   ,8,"Point Model"};
const Model gaussian_model = {eval_gaussian,9,"Gaussian Model"};
const Model gaussian_model_num_dev  = {eval_gaussian_num_dev,9,"Gaussian Model with Numerical Derivatives"};
const Model gaussian_model_testing = {eval_gaussian_testing,9,"Gaussian Model"};
