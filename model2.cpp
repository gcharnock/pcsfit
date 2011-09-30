
#include "model2.hpp"
#include "cuba.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <cmath>
#include <cassert>

using namespace std;

#define NDIM 3
#define NCOMP 1
#define EPSREL 1e-5
#define EPSABS 1000
#define VERBOSE 0
#define LAST 4
#define MINEVAL 100000
#define MAXEVAL 5000000

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
    }
    assert(false);
}



void eval_point(const double pm[8],double* value, double gradient[8]) {
    double x = pm[PARAM_X];
    double y = pm[PARAM_Y];
    double z = pm[PARAM_Z];

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
    double r7 = r5*r2;

    double inv12PiR5 = 1/(12*M_PI*r5);

    *value = inv12PiR5 *((r2-3*x2)*chi_1 + (z2-y2)*chi_2 + 6*(xy*chi_xy + xz*chi_xz + yz*chi_yz));
    gradient[PARAM_CHI1]  = (-2*x2+y2+z2)*inv12PiR5;
    gradient[PARAM_CHI2]  = (z2-y2)*inv12PiR5;
    gradient[PARAM_CHIXY] = xy*inv12PiR5;
    gradient[PARAM_CHIXZ] = xz*inv12PiR5;
    gradient[PARAM_CHIYZ] = yz*inv12PiR5;

    double A = (-2*x2+y2+z2)*chi_1 +(z2-y2)*chi_2 + 6*(xy*chi_xy + xz*chi_xz + yz*chi_yz);
    double fiveAOver12Pi7 = 5*A/(12*M_PI*r7);

    gradient[PARAM_X] = x*fiveAOver12Pi7 + (-4*x*chi_1          + 6*(y*chi_xy + z*chi_xz))*inv12PiR5;
    gradient[PARAM_Y] = y*fiveAOver12Pi7 + (2*y*(chi_1 - chi_2) + 6*(x*chi_xy + z*chi_yz))*inv12PiR5;
    gradient[PARAM_Z] = z*fiveAOver12Pi7 + (2*z*(chi_1 + chi_2) + 6*(x*chi_xz + y*chi_yz))*inv12PiR5;
}


struct Userdata {
    //For sending cuhre.
    const double pm[8];
    double stddev;

    //The position of 
    
    double xmax;
    double xmin;
               
    double ymax;
    double ymin;
               
    double zmax;
    double zmin;
};

//Since we need to pass this function to a C api it can't actually be
//a member function. We'll pass the this pointer explicitly to fake a
//member

int Integrand2(const int *ndim, const double xx[],
			  const int *ncomp, double ff[], void *voiduserdata) {
    Userdata* userdata = (Userdata*) voiduserdata;
    const double* pm = userdata->pm;
    double stddev = userdata->stddev;

	//Apply the transformation from the unit cube [0,1]^2. xp,yp,zp is
	//the dummy variable in molecular coordinates. The free variable
	//is intDet->xyz
	double xp = xx[0]*(userdata->xmax - userdata->xmin) + userdata->xmin;
	double yp = xx[1]*(userdata->ymax - userdata->ymin) + userdata->ymin;
	double zp = xx[2]*(userdata->zmax - userdata->zmin) + userdata->zmin;

	//Now evaulate rho at x-xp where x is the free variable
	double xp_m_x = xp - pm[PARAM_X];
	double yp_m_y = yp - pm[PARAM_Y];
	double zp_m_z = zp - pm[PARAM_Z];

    //Evauate the model
    double f;
    double gradient[8];
    eval_point(pm,&f,gradient);

    //We don't need the gradient of the gaussian function with respect
    //to the spaceial pramiters
    double a_coef = 1/(stddev*stddev);                                    assert(isfinite(a_coef));
    double theExp = exp(-a_coef*(xp_m_x*xp_m_x + yp_m_y*yp_m_y + zp_m_z*zp_m_z));
    double g =  pow(abs(a_coef)/M_PI,1.5)*theExp;  assert(isfinite(g));

    double dg = pow(M_PI,-1.5)*(1.5*pow(stddev,0.5)+2*pow(stddev,-1.5))*theExp;

	double result = g * f;
	//Scale the result
    double result_scaling = (userdata->xmax - userdata->xmin) *
		(userdata->ymax - userdata->ymin) *
		(userdata->zmax - userdata->zmin);
    
    //Pass the results and gradient
    ff[0] = result * result_scaling;
    for(unsigned long i=1;i<9;i++){ff[i] = gradient[i-1] * result_scaling;}
	ff[9] = dg * result_scaling;
            
    for(unsigned long i=0;i<10;i++){assert(isfinite(ff[i]));};

    return 0;
}

