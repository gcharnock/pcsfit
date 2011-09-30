
#include "model2.hpp"
#include "cuba.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <cmath>

using namespace std;

#define NDIM 3
#define NCOMP 1
#define EPSREL 1e-5
#define EPSABS 1000
#define VERBOSE 0
#define LAST 4
#define MINEVAL 100000
#define MAXEVAL 5000000

void eval_point(const PointModel* pm,double* value, PointModel* gradient) {
    double x = pm->metal.x;
    double y = pm->metal.y;
    double z = pm->metal.z;

    double chi_1 = pm->chi_1;
    double chi_2 = pm->chi_2;
    double chi_xy = pm->chi_xy;
    double chi_xz = pm->chi_xz;
    double chi_yz = pm->chi_yz;

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
    gradient->chi_1  = (-2*x2+y2+z2)*inv12PiR5;
    gradient->chi_2  = (z2-y2)*inv12PiR5;
    gradient->chi_xy = xy*inv12PiR5;
    gradient->chi_xz = xz*inv12PiR5;
    gradient->chi_yz = yz*inv12PiR5;

    double A = (-2*x2+y2+z2)*chi_1 +(z2-y2)*chi_2 + 6*(xy*chi_xy + xz*chi_xz + yz*chi_yz);
    double fiveAOver12Pi7 = 5*A/(12*M_PI*r7);

    gradient->metal.x =x*fiveAOver12Pi7 + (-4*x*chi_1          + 6*(y*chi_xy + z*chi_xz))*inv12PiR5;
    gradient->metal.y =y*fiveAOver12Pi7 + (2*y*(chi_1 - chi_2) + 6*(x*chi_xy + z*chi_yz))*inv12PiR5;
    gradient->metal.z =z*fiveAOver12Pi7 + (2*z*(chi_1 + chi_2) + 6*(x*chi_xz + y*chi_yz))*inv12PiR5;

}


struct Userdata {
    //For sending cuhre.
    const PointModel* pm;
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
    const PointModel* pm = userdata->pm;
    double stddev = userdata->stddev;

	//Apply the transformation from the unit cube [0,1]^2. xp,yp,zp is
	//the dummy variable in molecular coordinates. The free variable
	//is intDet->xyz
	double xp = xx[0]*(userdata->xmax - userdata->xmin) + userdata->xmin;
	double yp = xx[1]*(userdata->ymax - userdata->ymin) + userdata->ymin;
	double zp = xx[2]*(userdata->zmax - userdata->zmin) + userdata->zmin;

	//Now evaulate rho at x-xp where x is the free variable
	double xp_m_x = xp - pm->metal.x;
	double yp_m_y = yp - pm->metal.y;
	double zp_m_z = zp - pm->metal.z;

    //Evauate the model
    double f;
    PointModel gradient;
    eval_point(pm,&f,&gradient);
    

    //We don't need the gradient of the gaussian function, fortunatly
    double a_coef = 1/(stddev*stddev);                                    assert(isfinite(a_coef));
    double g = exp(-a_coef*(xp_m_x*xp_m_x+
							yp_m_y*yp_m_y+
							zp_m_z*zp_m_z)) * pow(abs(a_coef)/M_PI,1.5);  assert(isfinite(g));

	double result = g * f;
	//Scale the result
    double result_scaling = (userdata->xmax - userdata->xmin) *
		(userdata->ymax - userdata->ymin) *
		(userdata->zmax - userdata->zmin);
    
    //Pass the results and gradient
	ff[0] = result * result_scaling; //model

	ff[1] = result * result_scaling; //chi_1
	ff[2] = result * result_scaling; //chi_2
	ff[3] = result * result_scaling; //chi_xy
	ff[4] = result * result_scaling; //chi_xz
	ff[5] = result * result_scaling; //chi_yz

	ff[6] = result * result_scaling; //x
	ff[7] = result * result_scaling; //y
	ff[8] = result * result_scaling; //z

    assert(isfinite(ff[0]));

    return 0;
}

