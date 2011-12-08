
#include "model2.hpp"
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

double magic_polynomial(double t) {
    double t2 =  t*t;
    double t3 =  t2*t;
    double t4 =  t3*t;
    double t5 =  t4*t;
    double t6 =  t5*t;
    double t7 =  t6*t;
    double t8 =  t7*t;
    double t9 =  t8*t;
    double t10 = t9*t;
    
    return -t10*
        (-92378
         +755820*t
         -2771340* t2
         +5969040* t3
         -8314020* t4
         +7759752* t5
         -4849845* t6
         +1956240* t7
         -461890 * t8
         +48620  * t9);
}

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
    const static double EPSREL = 1e-6;
    const static double EPSABS = 0.00;
    const static int VERBOSE = 0;
    const static int LAST = 4;
    const static int MINEVAL = 1;
    const static int MAXEVAL = 5000000;
    const static int KEY = 7; //Use 11 point quadriture

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

void numerical_derivative(Vector3 evalAt,const Model* model,const double* params,double* gradient) {
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


void eval_point(Vector3 evalAt,const double* pm,double* value, double* gradient) {
    double x = evalAt.x - pm[PARAM_X];
    double y = evalAt.y - pm[PARAM_Y];
    double z = evalAt.z - pm[PARAM_Z];

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

        //I'm not sure yet what the meaning of the sign error. Added a - to the total expression to fix it.
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
	Vector3 evalAt;

    double reg_radius2;
};

double bump(double t)  {
    return abs(t) < 1 ? exp(-1/(1-t*t)) : 0;
}

//Since we need to pass this function to a C api it can't actually be
//a member function. We'll pass the this pointer explicitly to fake a
//member

int Integrand2(const double xx[],double ff[],int ncomp, void* void_userdata) {
    Userdata* userdata = (Userdata*)(void_userdata);

    const double* params = userdata->params;
    double stddev = abs(params[PARAM_STDDEV]);

    unsigned long point_size = userdata->point_model->size;

    double x = xx[0];
    double y = xx[1];
    double z = xx[2];
	double r2 = x*x + y*y + z*z;

    //A vector pointing from the centre of 1/r^3 to the center of the gaussian
    double singularity_x = userdata->evalAt.x - x - params[0];
    double singularity_y = userdata->evalAt.y - y - params[1];
    double singularity_z = userdata->evalAt.z - z - params[2];
    
    double singularity_r2 = singularity_x*singularity_x +
        singularity_y*singularity_y +
        singularity_z*singularity_z;

    //Evauate the model
    double f;
    double* gradient = ncomp == 1 ? NULL: (double*)alloca(point_size*sizeof(double));

    Vector3 x_minus_xprime(userdata->evalAt.x - x,
                           userdata->evalAt.y - y,
                           userdata->evalAt.z - z);

    userdata->point_model->modelf(x_minus_xprime,params,&f,gradient);

    //We don't need the gradient of the gaussian function with respect
    //to the spaceial pramiters
    double a_coef = 1/(stddev*stddev);                                    assert(isfinite(a_coef));

    //Potential optimisation: pull the calculation of this out of the integral
    double normalizer = pow(M_PI*stddev*stddev,-1.5);
    double theExp = exp(-a_coef*r2);
    double rho = normalizer*theExp;

    //Is any sort of regualisation needed?
    bool close      = singularity_r2 < 4*userdata->reg_radius2;
    bool very_close = singularity_r2 < userdata->reg_radius2;

    double another_r2;

    double another_x;
    double another_y;
    double another_z;

    double correction = 0;

    another_x  = (userdata->evalAt.x - params[0]);
    another_y  = (userdata->evalAt.y - params[1]);
    another_z  = (userdata->evalAt.z - params[2]);
    another_r2 = another_x*another_x + another_y*another_y + another_z*another_z;

    if(close) {
        if(very_close) {
            correction = 1;
        } else {
            double reg_radius = sqrt(userdata->reg_radius2);
            double singularity_r = sqrt(singularity_r2);

            double t =  (singularity_r - reg_radius)/(reg_radius);

            correction = magic_polynomial(t);
        }
    }
    double rho0 = normalizer*exp(-a_coef*(another_r2));

    assert(isfinite(rho));
    assert(isfinite(rho0));
    assert(isfinite(f));

    ff[0] = (rho-rho0*correction)*f;
    assert(isfinite(ff[0]));
    if(ncomp == 1) {
        //If we don't need a gradient, we can stop here.
        return 0;
    }

    double drhoBy_ds  = pow(M_PI,-1.5)*(2*pow(stddev,-6)*r2-3*pow(stddev,-4))*theExp;

    double theExp0 = exp(-a_coef*another_r2);
    double drhoBy_ds0 = pow(M_PI,-1.5)*(2*pow(stddev,-6)*another_r2-3*pow(stddev,-4))*theExp0;

    //double common_factor = -2*pow(M_PI,-3.0/2.0)/sqrt(stddev) * theExp0;
    //double drhoBydx_dx   = common_factor * singularity_x * (x-singularity_x);
    //double drhoBydy_dy   = common_factor * singularity_y * (y-singularity_y);
    //double drhoBydz_dz   = common_factor * singularity_z * (z-singularity_z);
        
    //double drho0 = drhoBydx_dx*another_x + drhoBydy_dy*another_y + drhoBydz_dz*another_z;

    for(unsigned long i=1;i<9;i++){
        ff[i] = (rho-rho0*correction) * gradient[i-1];
    }
    ff[9] = (drhoBy_ds-correction*drhoBy_ds0) * f;

    for(int i=0;i<ncomp;i++){assert(isfinite(ff[i]));}
    return 0;
}

//Is [c,d] a subset of [a,b]?
bool subset(double a,double b,double c,double d) {
    assert(a<b);
    assert(c<d);
    return (c<a && c<b) || (d<a && d<b);
}
bool intersect(double a,double b,double c,double d) {
    return subset(a,b,c,d) || subset(c,d,a,b);
}

void eval_gaussian(Vector3 evalAt,const double* params,double* value, double* gradient) {
    Userdata userdata;
    userdata.point_model = &point_model;
    userdata.params = params;
	userdata.evalAt = evalAt;

    double stddev = params[PARAM_STDDEV];

    IntegralBounds bounds;
    bounds.xmax =   12*stddev;
    bounds.xmin =  -12*stddev;
    
    bounds.ymax =   12*stddev;
    bounds.ymin =  -12*stddev;
    
    bounds.zmax =   12*stddev;
    bounds.zmin =  -12*stddev;

    //Make sure that the region where regualization is applied is
    //wholey included or wholey excluded.

    double reg_radius    = stddev/3.0;
    userdata.reg_radius2 = reg_radius*reg_radius;
    unsigned long ncomp = gradient == NULL ? 1 : 10;

	double* integral = (double*)alloca(ncomp*sizeof(double));

	cuhreIntegrate(Integrand2,&bounds,ncomp,integral,(void*)&userdata);

    *value = integral[0];

    if(gradient != NULL) {
        memcpy(gradient,integral+1,10*sizeof(double));
    }
}


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
		dataset->nuclei[i] = Vector3(x + 10*rand(prng),y + 10*rand(prng),z + 10*rand(prng));
        model.modelf(dataset->nuclei[i],params,&(dataset->vals[i]),NULL);
    }
}

void eval_gaussian_testing(Vector3 evalAt,const double* params,double* value, double* gradient) {
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

                Vector3 x_minus_xprime(evalAt.x - x,
                                       evalAt.y - y,
                                       evalAt.z - z);
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

    Integrand2(xx,ff,1,(void*)&userdata_copy);

    for(unsigned long i = 0; i < size; i++) {
        double param = userdata->params[i];

		double h     = abs(param*0.00001);
		double result_plus  = numeric_limits<double>::quiet_NaN();
        double result_minus = numeric_limits<double>::quiet_NaN();

		params_mutable[i] = param + h;
        Integrand2(xx,&result_plus,1,(void*)&userdata_copy);

		params_mutable[i] = param - h;
        Integrand2(xx,&result_minus,1,(void*)&userdata_copy);

		params_mutable[i] = param;

		ff[i+1] = (result_plus-result_minus)/(2*h);
        //if(i == 8 && abs(xx[2]) < 0.01) cout << xx[0] << " " << xx[1] << " " <<  ff[i+1] << endl;
	}

    return 0;
}


void eval_gaussian_num_dev(Vector3 evalAt,const double* params,double* value, double* gradient) {
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
        cuhreIntegrate(Integrand2     ,&bounds,ncomp,integral,(void*)&userdata);
    } else {
        cuhreIntegrate(IntegrandNumDev,&bounds,ncomp,integral,(void*)&userdata);
    }

    *value = integral[0];

    if(gradient != NULL) {
        memcpy(gradient,integral+1,9*sizeof(double));
    }
}

const Model point_model    = {eval_point   ,8,"Point Model"};
const Model gaussian_model = {eval_gaussian,9,"Gaussian Model"};
const Model gaussian_model_num_dev  = {eval_gaussian_num_dev,9,"Gaussian Model with Numerical Derivatives"};
const Model gaussian_model_testing = {eval_gaussian_testing,9,"Gaussian Model"};
