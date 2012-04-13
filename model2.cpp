
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


void cuhreIntegrate(IntegrandF f,
                    IntegralBounds* bounds,
                    unsigned long ncomp,
                    double* integral,
                    void* data,
                    double epsRel,
                    double epsAbs) {
    const static int NDIM = 3;
    const static int VERBOSE = 0;
    const static int LAST = 4;
    const static int MINEVAL = 1;
    const static int MAXEVAL = 200000000;
    const static int KEY = 11; //Use 11 point quadriture

	int nregions, neval, fail;
	double* error = (double*)alloca(ncomp * sizeof(double));
    double* prob  = (double*)alloca(ncomp * sizeof(double));

	bounds->integrand=f;
    bounds->data=data;

	Cuhre(NDIM, ncomp, cuhreIntegrand, (void*)bounds,
		  epsRel, epsAbs, VERBOSE | LAST,
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

void numerical_derivative(Vec3d evalAt,
                          const Model* model,
                          const double* params,
                          double* gradient,
                          const ModelOptions* modelOptions) {
	double* params_mutable = (double*)alloca(model->size*sizeof(double));
	memcpy(params_mutable,params,model->size*sizeof(double));

    for(unsigned long i = 0;i<model->size;i++) {
		double h     = abs(params[i]*0.001);
		double result_plus,result_minus;


		params_mutable[i] = params[i] + h;
		model->modelf(evalAt,params_mutable,&result_plus ,NULL,modelOptions);

		params_mutable[i] = params[i] - h;
		model->modelf(evalAt,params_mutable,&result_minus,NULL,modelOptions);

		params_mutable[i] = params[i];

		gradient[i] = (result_plus-result_minus)/(2*h);
	}
}


void eval_point(Vec3d evalAt,const double* pm,double* value, double* gradient,const ModelOptions*) {
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


    //The chi tensor has units of m^3 x 10^(-32), but we need it in
    //simple m^3 so divide by 10^-32 but our distances are in
    //angrstroms and shifts are in ppm, chi it needs to be 10^36
    //bigger to compensate, hence the net factor of 10^4
    double chi_1 =  pm[PARAM_CHI1]  * 1e4; 
    double chi_2 =  pm[PARAM_CHI2]  * 1e4; 
    double chi_xy = pm[PARAM_CHIXY] * 1e4;
    double chi_xz = pm[PARAM_CHIXZ] * 1e4;
    double chi_yz = pm[PARAM_CHIYZ] * 1e4;

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

    double A = (3*z2-r2)*chi_1 + (x2-y2)*chi_2 + 6*(xy*chi_xy + xz*chi_xz + yz*chi_yz);
    *value = inv12PiR5 * A; assert(isfinite(*value));

    if(gradient != NULL) {
        gradient[PARAM_CHI1]  = inv12PiR5*(3*z2-r2)    / 1e-4;
        gradient[PARAM_CHI2]  = inv12PiR5*(x2-y2)      / 1e-4;
        gradient[PARAM_CHIXY] = 6*inv12PiR5*xy         / 1e-4;
        gradient[PARAM_CHIXZ] = 6*inv12PiR5*xz         / 1e-4;
        gradient[PARAM_CHIYZ] = 6*inv12PiR5*yz         / 1e-4;


        double fiveAOver12Pi7 = 5*A/(12*M_PI*r5*r2);

        //There is a minus here because this is a deriative with
        //respect to the position of the metal, which appears as
        //sigma(r - r_m) where sigma is a model assuming the metal
        //sits at the origin and r_m is the location of the metal
        gradient[PARAM_X] = -(-x*fiveAOver12Pi7 + (-2*x*(chi_1 - chi_2) + 6*(y*chi_xy + z*chi_xz))*inv12PiR5);
        gradient[PARAM_Y] = -(-y*fiveAOver12Pi7 + (-2*y*(chi_1 + chi_2) + 6*(x*chi_xy + z*chi_yz))*inv12PiR5);
        gradient[PARAM_Z] = -(-z*fiveAOver12Pi7 + ( 4*z*chi_1           + 6*(x*chi_xz + y*chi_yz))*inv12PiR5);
	
        for(unsigned long i = 0;i<8;i++) {
            assert(isfinite(gradient[i]));
            assert(isfinite(gradient[i]*gradient[i]));
        }
    }
}

enum INTEGRAL_TYPE {
    COMPUTE_ANALYTIC,
    COMPUTE_X,
    COMPUTE_Y,
    COMPUTE_Z,
    COMPUTE_STDDEV
};

struct Userdata {
    //For sending cuhre.
    const Model* point_model;
    const double* params; 
	Vec3d evalAt;

    //Fintie differencing step size
    double stddev_h; 
    double position_h;

    double reg_radius2;

    INTEGRAL_TYPE integral_type;

    const ModelOptions* modelOptions;
};




int parse_params_file(const std::string& filename,const Model** model,std::vector<double>* params) {
    ifstream fin(filename.c_str());
    if(!fin.is_open()) {
        return PARAM_FILE_NOT_FOUND;
    }
    string model_name;
    fin >> model_name;
    if(model_name == "point" || model_name == "euler_point") {
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

    if(model_name == "euler_point") {
        //The params are actual x,y,z,ax, rh, alpha, beta, gamma so we
        //need to convert the last three

        AxRhomTensor axRhomTensor(params->at(3),params->at(4),params->at(5)/180*M_PI,params->at(6)/180*M_PI,params->at(7)/180*M_PI);
        Tensor converted = axRhomToTensor(axRhomTensor);
        (*params)[3] = converted.chi_1;
        (*params)[4] = converted.chi_2;
        (*params)[5] = converted.chi_xy;
        (*params)[6] = converted.chi_xz;
        (*params)[7] = converted.chi_yz;
    }

    return PARSE_SUCESS;
}


void random_data(PRNG& prng,
                 const Model& model,
                 const double* params,
                 unsigned long natoms,
                 Dataset* dataset,
                 const ModelOptions* modelOptions) {
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
        model.modelf(dataset->nuclei[i],params,&(dataset->vals[i]),NULL,modelOptions);
    }
}

//================================================================================//

double integrand_kernel(Userdata* userdata,
                        Vec3d xyz,
                        const double* params,
                        double* gradient) {

    for(ulong i = 0; i < POINT_SIZE+1; i++) {
        assert(isfinite(params[i]));
    }

    Vec3d metal = Vec3d(params);

    //The location of the singularity
    Vec3d singularity = userdata->evalAt - metal;

    assert(isfinite(xyz[0]));
    assert(isfinite(xyz[1]));
    assert(isfinite(xyz[2]));

    assert(isfinite(singularity[0]));
    assert(isfinite(singularity[1]));
    assert(isfinite(singularity[2]));

    double singularity_r2 = singularity.r2();
    double stddev = abs(params[PARAM_STDDEV]);
    double stddev2 = stddev*stddev;


    //Potential optimisation: pull the calculation of these out of the integral
    double a_coef = 1/(stddev*stddev);   assert(isfinite(a_coef));
    double normalizer = pow(M_PI*stddev*stddev,-1.5);
    double theExp0 = exp(-a_coef*(singularity_r2));

    //rho0 and its first derivative
    double rho0 = normalizer*theExp0;
    Vec3d rho0_grad = singularity *(-2/stddev2 * rho0);

    //Read the dummy variable from xx[]
	double r2 = xyz.r2();

    //A vector pointing from the centre of 1/r^3 to the center of the gaussian
    Vec3d expand_around = xyz - singularity;
    double expand_r2 = expand_around.r2();


    //Decide if any sort of regualisation needed?
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



    //Evauate the model
    double f;
    double* point_gradient = gradient == NULL ? NULL : (double*)alloca(POINT_SIZE*sizeof(double));

    Vec3d x_minus_xprime = userdata->evalAt - xyz;

    userdata->point_model->modelf(x_minus_xprime,params,&f,point_gradient,userdata->modelOptions);

    double theExp = exp(-a_coef*r2);
    double rho = normalizer*theExp;

    double toSub =  correction*(rho0 + rho0_grad.dot(expand_around));

    if(abs((rho-toSub)/rho) < 1e-10) {
        //cout << " Warning, possible loss of precision, rho = "
        //<< rho << " toSub = " << toSub << "  diff = " << (rho-toSub) << endl;

        //cout << " computed diff = " << normalizer*gaussian_error_term_one(r,a,stddev) << endl;
    }

    if(gradient != NULL) {
        //Ignore the spacial gradients, because they don't behave well
        //under integration.
        for(unsigned long i=0;i<5;i++){
            gradient[i] = (rho-toSub) * point_gradient[i+3];
        }
    }

    double retVal = (rho - toSub*correction)*f;
    assert(isfinite(retVal));
    return retVal;
}


int Integrand(const double xx[],double ff[],int ncomp, void* void_userdata) {
    /*
      The meaning of the output vector depends on the value of
      userdata->integral_type
      
      COMPUTE_ANALYTIC
      ff[0]     - The function value
      ff[1-5]   - Tensor component gradients

      COMPUTE_X
      ff[0-3] - x centeral differences, running from f(s-h) to f(s+h)

      COMPUTE_Y
      ff[0-3] - y centeral differences, running from f(s-h) to f(s+h)

      COMPUTE_Z
      ff[0-3] - z centeral differences, running from f(s-h) to f(s+h)

      COMPUTE_STDDEV
      ff[0-7] - s centeral differences, running from f(s-4*h) to f(s+4*h)
    */

    Userdata* userdata = (Userdata*)(void_userdata);

    const double* params = userdata->params;

    Vec3d xyz = Vec3d(xx);

    //Evauate the model
    double* gradient = ncomp == 1 ? NULL: (double*)alloca(POINT_SIZE*sizeof(double));

    if(userdata->integral_type == COMPUTE_ANALYTIC) {
        ff[0] = integrand_kernel(userdata,xyz,params,gradient);
        assert(isfinite(ff[0]));
    }

    if(ncomp == 1) {
        assert(userdata->integral_type == COMPUTE_ANALYTIC);
        //If we don't need a gradient, we can stop here.
        return 0;
    }
    if(userdata->integral_type == COMPUTE_ANALYTIC) {
        assert(ncomp = 6);
        //Copy the analytic gradients over
        memcpy(ff+1,gradient,5*sizeof(double));
    }
    
    //The stddev central differences
    double* mute_params = (double*)alloca((POINT_SIZE+1)*sizeof(double));
    memcpy(mute_params,params,(POINT_SIZE+1)*sizeof(double));

    if(userdata->integral_type == COMPUTE_STDDEV) {

        mute_params[PARAM_STDDEV] = params[PARAM_STDDEV] - userdata->stddev_h * 4;
        ff[0] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-4h
        mute_params[PARAM_STDDEV] = params[PARAM_STDDEV] - userdata->stddev_h * 3;
        ff[1] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-3h
        mute_params[PARAM_STDDEV] = params[PARAM_STDDEV] - userdata->stddev_h * 2;
        ff[2] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-2h
        mute_params[PARAM_STDDEV] = params[PARAM_STDDEV] - userdata->stddev_h * 1;
        ff[3] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-1h

        mute_params[PARAM_STDDEV] = params[PARAM_STDDEV] + userdata->stddev_h * 1;
        ff[4] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 1h
        mute_params[PARAM_STDDEV] = params[PARAM_STDDEV] + userdata->stddev_h * 2;
        ff[5] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 2h
        mute_params[PARAM_STDDEV] = params[PARAM_STDDEV] + userdata->stddev_h * 3;
        ff[6] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 3h
        mute_params[PARAM_STDDEV] = params[PARAM_STDDEV] + userdata->stddev_h * 4;
        ff[7] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 4h

        mute_params[PARAM_STDDEV] = params[PARAM_STDDEV];

    } else if(userdata->integral_type == COMPUTE_X) {

        //Vary x
        mute_params[PARAM_X] = params[PARAM_X] - userdata->position_h*3;
        ff[0] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-3h
        mute_params[PARAM_X] = params[PARAM_X] - userdata->position_h*2;
        ff[1] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-2h
        mute_params[PARAM_X] = params[PARAM_X] - userdata->position_h;
        ff[2] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-1h

        mute_params[PARAM_X] = params[PARAM_X] + userdata->position_h;
        ff[3] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 1h
        mute_params[PARAM_X] = params[PARAM_X] + userdata->position_h*2;
        ff[4] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 2h
        mute_params[PARAM_X] = params[PARAM_X] + userdata->position_h*3;
        ff[5] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 3h
    
        mute_params[PARAM_X] = params[PARAM_X];

    } else if(userdata->integral_type == COMPUTE_Y) {


        //Vary y
        mute_params[PARAM_Y] = params[PARAM_Y] - userdata->position_h*3;
        ff[0] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-3h
        mute_params[PARAM_Y] = params[PARAM_Y] - userdata->position_h*2;
        ff[1] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-2h
        mute_params[PARAM_Y] = params[PARAM_Y] - userdata->position_h;
        ff[2] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-1h

        mute_params[PARAM_Y] = params[PARAM_Y] + userdata->position_h;
        ff[3] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 1h
        mute_params[PARAM_Y] = params[PARAM_Y] + userdata->position_h*2;
        ff[4] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 2h
        mute_params[PARAM_Y] = params[PARAM_Y] + userdata->position_h*2;
        ff[5] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 3h

        mute_params[PARAM_Y] = params[PARAM_Y];

    } else if(userdata->integral_type == COMPUTE_Z) {

        //Vary z
        mute_params[PARAM_Z] = params[PARAM_Z] - userdata->position_h*3;
        ff[0] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-3h
        mute_params[PARAM_Z] = params[PARAM_Z] - userdata->position_h*2;
        ff[1] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-2h
        mute_params[PARAM_Z] = params[PARAM_Z] - userdata->position_h;
        ff[2] = integrand_kernel(userdata,xyz,mute_params,gradient);  //-1h

        mute_params[PARAM_Z] = params[PARAM_Z] + userdata->position_h;
        ff[3] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 1h
        mute_params[PARAM_Z] = params[PARAM_Z] + userdata->position_h*2;
        ff[4] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 2h
        mute_params[PARAM_Z] = params[PARAM_Z] + userdata->position_h*3;
        ff[5] = integrand_kernel(userdata,xyz,mute_params,gradient);  // 3h

        mute_params[PARAM_Z] = params[PARAM_Z];
    }

    for(int i=0;i<ncomp;i++){
        assert(isfinite(ff[i]));
    }
    return 0;
}


void eval_gaussian(Vec3d evalAt,const double* params,double* value, double* gradient,const ModelOptions* modelOptions) {
    Userdata userdata;
    userdata.point_model = &point_model;
    userdata.params = params;
	userdata.evalAt = evalAt;
    userdata.modelOptions = modelOptions;

    double relError = modelOptions->relError;
    double absError = modelOptions->absError;

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

    //The size of the integral result buffer to allocate. If we don't
    //care about gradients, we will only be calling cuhreIntegrate
    //once for a single integral, so set this to 1. Otherwise the
    //maximum value is 8, which is used by the 4th order finite
    //difference scheme.
    unsigned long max_ncomp = gradient == NULL ? 1 : 8;

	double* integral = (double*)alloca(max_ncomp*sizeof(double));

    userdata.stddev_h   = 0.002;
    userdata.position_h = 0.002;

    userdata.integral_type = COMPUTE_ANALYTIC;
	cuhreIntegrate(Integrand,
                   &bounds,
                   gradient == NULL ? 1 : 6, //6 is for value + 5 tensor gradients
                   integral,
                   (void*)&userdata,
                   relError,
                   absError);

    *value = integral[0];

    if(gradient != NULL) {
        //Copy the gradients we calculate analytically before reusing
        //the integral buffer
        gradient[PARAM_CHI1 ] = integral[1];
        gradient[PARAM_CHI2 ] = integral[2];
        gradient[PARAM_CHIXY] = integral[3];
        gradient[PARAM_CHIXZ] = integral[4];
        gradient[PARAM_CHIYZ] = integral[5];


        //Calculate and copy the calculated finite difference values for x
        userdata.integral_type = COMPUTE_X;
        cuhreIntegrate(Integrand,&bounds,6,integral,(void*)&userdata,relError,absError);

        double val_x3m = integral[0];
        double val_x2m = integral[1];
        double val_xm  = integral[2];
        double val_xp  = integral[3];
        double val_x2p = integral[4];
        double val_x3p = integral[5];

        //Calculate and copy the calculated finite difference values for y
        userdata.integral_type = COMPUTE_Y;
        cuhreIntegrate(Integrand,&bounds,6,integral,(void*)&userdata,relError,absError);

        double val_y3m = integral[0];
        double val_y2m = integral[1];
        double val_ym  = integral[2];
        double val_yp  = integral[3];
        double val_y2p = integral[4];
        double val_y3p = integral[5];

        //Calculate and copy the calculated finite difference values for z
        userdata.integral_type = COMPUTE_Z;
        cuhreIntegrate(Integrand,&bounds,6,integral,(void*)&userdata,relError,absError);

        double val_z3m = integral[0];
        double val_z2m = integral[1];
        double val_zm  = integral[2];
        double val_zp  = integral[3];
        double val_z2p = integral[4];
        double val_z3p = integral[5];

        //Calculate and copy the calculated finite difference values for stddev
        userdata.integral_type = COMPUTE_STDDEV;
        cuhreIntegrate(Integrand,&bounds,8,integral,(void*)&userdata,relError,absError);

        double val_4m = integral[0];
        double val_3m = integral[1];
        double val_2m = integral[2];
        double val_m  = integral[3];

        double val_p  = integral[4];
        double val_2p = integral[5];
        double val_3p = integral[6];
        double val_4p = integral[7];

        
        //double fd1gradient = (val_p-val_m)/(2*h);
        //double fd2gradient = (val_2m - 8*val_m + 8*val_p - val_2p)/(12*h);
        //double fd3gradient = (-val_3m + 9*val_2m - 45*val_m + 45*val_p - 9*val_2p + val_3p)/(60*h);

        double h = userdata.position_h;

        //double fd1xgradient = (val_xp-val_xm)/(2*h);
        //double fd2xgradient = (val_x2m - 8*val_xm + 8*val_xp - val_x2p)/(12*h);

        gradient[PARAM_X] = (-val_x3m + 9*val_x2m - 45*val_xm + 45*val_xp - 9*val_x2p + val_x3p)/(60*h);
        gradient[PARAM_Y] = (-val_y3m + 9*val_y2m - 45*val_ym + 45*val_yp - 9*val_y2p + val_y3p)/(60*h);
        gradient[PARAM_Z] = (-val_z3m + 9*val_z2m - 45*val_zm + 45*val_zp - 9*val_z2p + val_z3p)/(60*h);

        /*cout << "fd1xgradient = " << fd1xgradient 
             << " fd2xgradient = " << fd2xgradient
             << " fd3xgradient = " << gradient[PARAM_X] << endl;*/

        h = userdata.stddev_h;
        gradient[PARAM_STDDEV] = ( val_4m/280 - 4*val_3m/105 + val_2m/5 - 4*val_m/5
                                   -val_4p/280 + 4*val_3p/105 - val_2p/5 + 4*val_p/5)/h;

    }
}

//================================================================================//

void eval_gaussian_testing(Vec3d evalAt,
                           const double* params,
                           double* value,
                           double* gradient,
                           const ModelOptions* modelOptions) {
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
                eval_point(x_minus_xprime,params,&f,NULL,modelOptions);

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


void eval_gaussian_num_dev(Vec3d evalAt,const double* params,double* value, double* gradient,const ModelOptions* modelOptions) {
    Userdata userdata;
    userdata.point_model = &point_model;
    userdata.params = params;
	userdata.evalAt = evalAt;
    userdata.modelOptions = modelOptions;
    double stddev = params[PARAM_STDDEV];

    double relError = modelOptions->relError;
    double absError = modelOptions->absError;

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
        cuhreIntegrate(Integrand      ,&bounds,ncomp,integral,(void*)&userdata,relError,absError);
    } else {
        cuhreIntegrate(IntegrandNumDev,&bounds,ncomp,integral,(void*)&userdata,relError,absError);
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
const Model gaussian_model_testing =  {eval_gaussian_testing,9,"Gaussian Model"};
