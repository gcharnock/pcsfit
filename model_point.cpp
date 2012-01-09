

#include "model_common.hpp"
#include <cmath>

using namespace std;

const Model point_model    = {eval_point   ,8,"Point Model"};

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
