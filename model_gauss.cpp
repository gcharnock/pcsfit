
#include "model_common.hpp"


#include "maths.hpp"
#include "integrate.hpp"
#include "pointdev.hpp"

#include <cmath>

using namespace std;


int Integrand(const double xx[],double ff[],int ncomp, void* void_userdata) {
    Userdata* userdata = (Userdata*)(void_userdata);

    const double* params = userdata->params;
    double stddev = abs(params[PARAM_STDDEV]);

    unsigned long point_size = userdata->point_model->size;

    double x = xx[0];
    double y = xx[1];
    double z = xx[2];
	double r2 = x*x + y*y + z*z;

    //A vector pointing from the centre of 1/r^3 to the center of the gaussian
    double singularity_x = userdata->evalAt.x() - x - params[0];
    double singularity_y = userdata->evalAt.y() - y - params[1];
    double singularity_z = userdata->evalAt.z() - z - params[2];
    
    double singularity_r2 = singularity_x*singularity_x +
        singularity_y*singularity_y +
        singularity_z*singularity_z;

    //Evauate the model
    double f;
    double* gradient = ncomp == 1 ? NULL: (double*)alloca(point_size*sizeof(double));

    Vec3d x_minus_xprime(userdata->evalAt.x() - x,
                         userdata->evalAt.y() - y,
                         userdata->evalAt.z() - z);

    userdata->point_model->modelf(x_minus_xprime,params,&f,gradient);

    //We don't need the gradient of the gaussian function with respect
    //to the spaceial pramiters
    double a_coef = 1/(stddev*stddev);                       assert(isfinite(a_coef));

    //Potential optimisation: pull the calculation of this out of the integral
    double normalizer = pow(M_PI*stddev*stddev,-1.5);
    double theExp = exp(-a_coef*r2);
    double rho = normalizer*theExp;

    //Is any sort of regualisation needed?
    bool close      = singularity_r2 < 4*userdata->reg_radius2;
    bool very_close = singularity_r2 <   userdata->reg_radius2;

    double another_r2;

    double another_x;
    double another_y;
    double another_z;

    double correction = 0;

    another_x  = (userdata->evalAt.x() - params[0]);
    another_y  = (userdata->evalAt.y() - params[1]);
    another_z  = (userdata->evalAt.z() - params[2]);
    another_r2 = another_x*another_x + another_y*another_y + another_z*another_z;

    if(close) {
        if(very_close) {
            correction = 1;
        } else {
            double reg_radius = sqrt(userdata->reg_radius2);
            double singularity_r = sqrt(singularity_r2);

            double t =  (singularity_r - reg_radius)/(reg_radius);

            correction = magic_polynomial(1-t);
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

void eval_gaussian(Vec3d evalAt,const double* params,double* value, double* gradient) {
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

	cuhreIntegrate(Integrand,&bounds,ncomp,integral,(void*)&userdata);

    *value = integral[0];

    if(gradient != NULL) {
        memcpy(gradient,integral+1,10*sizeof(double));
    }
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
        cuhreIntegrate(Integrand     ,&bounds,ncomp,integral,(void*)&userdata);
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
const Model gaussian_model_num_dev  = {eval_gaussian_num_dev,9,"Gaussian Model with Numerical Derivatives"};

