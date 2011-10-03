
#include "model2.hpp"
#include "cuba.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <cstring>
#include <cmath>
#include <cassert>
#include "panic.hpp"

using namespace std;

struct IntegralBounds;

typedef int (*IntegrandF) (const double xx[],double ff[], IntegralBounds *bounds);

struct IntegralBounds {
    double xmax;
    double xmin;
               
    double ymax;
    double ymin;
               
    double zmax;
    double zmin;

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

	bounds->integrand(xxprime,ff,bounds);

    for(unsigned long i=0;i<10;i++){assert(isfinite(ff[i]));}

    return 0;
}

void cuhreIntegrate(IntegrandF f,IntegralBounds* bounds,double* integral) {
    const static int NDIM = 3;
    const static int NCOMP = 10;
    const static double EPSREL = 1e-2;
    const static double EPSABS = 0;
    const static int VERBOSE = 0;
    const static int LAST = 4;
    const static int MINEVAL = 1;
    const static int MAXEVAL = 5000000;

    const static int KEY = 11;

	int nregions, neval, fail;
	double error[10], prob[10];

	bounds->integrand=f;

	Cuhre(NDIM, NCOMP, cuhreIntegrand, (void*)bounds,
		  EPSREL, EPSABS, VERBOSE | LAST,
		  MINEVAL, MAXEVAL, KEY,
		  &nregions, &neval, &fail, integral, error, prob);

    double result_scaling = (bounds->xmax - bounds->xmin) *
		(bounds->ymax - bounds->ymin) *
		(bounds->zmax - bounds->zmin);
    for(unsigned long i = 0; i < 10; i++){integral[i]*=result_scaling;}

    for(unsigned long i = 0; i < 10;i++) {
        assert(isfinite(integral[i]));
    }
    if(fail != 0)  {
        cerr << "Warning, an integral failed to converge" << endl;
    }
}


std::string name_param(POINT_PARAM param) {
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
    PANIC();
    return "";
}

void numerical_derivative(Vector3 evalAt,double* model,ModelF modelf,unsigned long nparams,double * gradient) {
    assert(nparams < MAX_PARAMS);
    double fake_gradient[MAX_PARAMS];

    for(unsigned long i = 0;i<8;i++) {
		double h     = abs(model[i]*0.000001);
		double result_plus,result_minus;

		double original = model[i];

		model[i] = original + h;
		modelf(evalAt,model,&result_plus ,fake_gradient);
		model[i] = original - h;
		modelf(evalAt,model,&result_minus,fake_gradient);

		model[i] = original;

		gradient[i] = (result_plus-result_minus)/(2*h);
	}
}


void eval_point(Vector3 evalAt,double pm[8],double* value, double gradient[8]) {
    double x = evalAt.x - pm[PARAM_X];
    double y = evalAt.y - pm[PARAM_Y];
    double z = evalAt.z - pm[PARAM_Z];

	if(x == 0 && y == 0 && z == 0) {
		*value = 0;
		for(unsigned long i = 0;i < 8; i++){gradient[i] = 0;}
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
    *value = inv12PiR5 * A;

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
	
	for(unsigned long i = 0;i<8;i++) {assert(isfinite(gradient[i]));}
}

void eval_point_ND(Vector3 evalAt,double pm[8],double* value) {
    double x = evalAt.x - pm[PARAM_X];
    double y = evalAt.y - pm[PARAM_Y];
    double z = evalAt.z - pm[PARAM_Z];

	if(x == 0 && y == 0 && z == 0) {
		*value = 0;
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
    *value = inv12PiR5 * A;
	
	assert(isfinite(*value));
}



struct Userdata : public IntegralBounds {
    //For sending cuhre.
    double* pm; //Should be of length 8
    double  stddev;
	Vector3 evalAt;
};

//Since we need to pass this function to a C api it can't actually be
//a member function. We'll pass the this pointer explicitly to fake a
//member

int Integrand2(const double xx[],double ff[], IntegralBounds* bounds) {
    Userdata* userdata = static_cast<Userdata*>(bounds);

    double* pm = userdata->pm;
    double stddev = userdata->stddev;

	if(pm[PARAM_X] == 0.0 && pm[PARAM_Y] == 0.0 && pm[PARAM_Z] == 0.0) {
		for(unsigned long i = 0; i < 10; i++){ff[i] = 0.0;}
		return 0;
	}

    //Evauate the model
    double f;
    double gradient[8];
    eval_point(userdata->evalAt,pm,&f,gradient);

    //We don't need the gradient of the gaussian function with respect
    //to the spaceial pramiters
    double a_coef = 1/(stddev*stddev);                                    assert(isfinite(a_coef));
    double theExp = exp(-a_coef*(xx[0]*xx[0] + xx[1]*xx[1] + xx[2]*xx[2]));
    double rho =  theExp;  assert(isfinite(rho));

    double drho = pow(M_PI,-1.5)*(1.5*pow(stddev,0.5)+2*pow(stddev,-1.5))*theExp;

    //Pass the results and gradient
    ff[0] = rho * f;
    for(unsigned long i=1;i<9;i++){ff[i] = rho * gradient[i-1];}
	ff[9] = drho * f;

    return 0;
}

void eval_gaussian(Vector3 evalAt,double* model,double* value, double gradient[9]) {
    Userdata userdata;
    userdata.pm = model;
    userdata.stddev = model[PARAM_STDDEV];
	userdata.evalAt = evalAt;

    userdata.xmax =  5*userdata.stddev;
	userdata.xmin = -5*userdata.stddev;
    
    userdata.ymax =  5*userdata.stddev;
    userdata.ymin = -5*userdata.stddev;
    
    userdata.zmax =  5*userdata.stddev;
    userdata.zmin = -5*userdata.stddev;

	double integral[10];

	cuhreIntegrate(Integrand2,static_cast<IntegralBounds*>(&userdata),integral);

	//Scale the result
    double a_coef = 1/(userdata.stddev*userdata.stddev);         assert(isfinite(a_coef));
    double normalizer = pow(a_coef/M_PI,1.5);

    for(unsigned long i = 0; i < 10; i++){integral[i]*=normalizer;}


    *value = integral[0];
    memcpy(gradient,integral+1,9*sizeof(double));
}
