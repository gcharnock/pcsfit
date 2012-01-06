
#include "pointdev.hpp"
#include "model2.hpp"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "maths.hpp"

using namespace std;
using namespace boost;





//Evauates the numerator and derivatives of the point model
/*
  Derivatives come out in this order
  {0,0,0},
  
  {1,0,0},
  {0,1,0},
  {0,0,1},

  {1,1,0},
  {1,0,1},
  {0,1,1},

  {2,0,0},
  {0,2,0},
  {0,0,2}
*/

void eval_point_numerator(Vector3 evalAt,const double* pm,double out[9]) {
    double x = evalAt.x - pm[PARAM_X];
    double y = evalAt.y - pm[PARAM_Y];
    double z = evalAt.z - pm[PARAM_Z];

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

    out[0] = (r2-3*x2)*chi_1 + (z2-y2)*chi_2 + 6*(xy*chi_xy + xz*chi_xz + yz*chi_yz);

    out[1] = -4*x*chi_1 +             6*(y*chi_xy + z*chi_xz);
    out[2] =  2*y*chi_1 - 2*y*chi_2 + 6*(x*chi_xy + z*chi_yz);
    out[3] =  2*z*chi_1 + 2*z*chi_2 + 6*(x*chi_xz + y*chi_yz);

    out[4] = 6*chi_xy;
    out[5] = 6*chi_xz;
    out[6] = 6*chi_yz;

    out[7] = -4*chi_1;   
    out[8] =  2*chi_1 - 2*chi_2;
    out[9] =  2*chi_1 + 2*chi_2;
    cout << "Azz = " << out[9] << endl;

    for(unsigned long i = 0;i<10;i++) {
        assert(isfinite(out[i]));
    }
}



double u_nmp_over_v(ulong n,ulong m,ulong p, double x,double y,double z) {
    MNL u = MNL::one();
    for(ulong i = 0; i<n; i++) {
        u = u_plus_one_recusion(u,0,i);
    }
    for(ulong i = 0; i<m; i++) {
        u = u_plus_one_recusion(u,1,i+n);
    }
    for(ulong i = 0; i<p; i++) {
        u = u_plus_one_recusion(u,2,i+n+m);
    }

    ulong nTot = n + m + p;

    double vn = pow(x*x+y*y+z*z, (2.0*nTot+5.0)/2.0);

    double xyz[3];
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
    
    cout << u << endl;

    return u.eval(xyz)/vn;
}

double eval_point_model_dev_xyz(const double* params,ulong ndx, ulong ndy, ulong ndz,Vector3 evalAt) {
    double x = evalAt.x - params[PARAM_X];
    double y = evalAt.y - params[PARAM_Y];
    double z = evalAt.z - params[PARAM_Z];

    double numerator_devs[10];
    eval_point_numerator(evalAt,params,numerator_devs);

    //There are ten possible partitions of 2, so we will just write
    //them out explicitly

    static const ulong two_partitions[10][3] = {
        {0,0,0},

        {1,0,0},
        {0,1,0},
        {0,0,1},

        {1,1,0},
        {1,0,1},
        {0,1,1},

        {2,0,0},
        {0,2,0},
        {0,0,2}
    };

    //Each of these partitions has an associated multinomial
    //coefficent NB: This must be in the same order as the
    //two_partitions array
    ulong multinomial_coef[10];

    multinomial_coef[0] = 1;

    multinomial_coef[1] = ndx;
    multinomial_coef[2] = ndy;
    multinomial_coef[3] = ndz;

    multinomial_coef[4] = ndx*ndy;
    multinomial_coef[5] = ndx*ndz;
    multinomial_coef[6] = ndy*ndz;

    multinomial_coef[7] = ndx*(ndx-1)/2;
    multinomial_coef[8] = ndy*(ndy-1)/2;
    multinomial_coef[9] = ndz*(ndz-1)/2;
    
    //Compute the derivative
    double total = 0;

    for(ulong i = 0;i<10;i++) {
        if(multinomial_coef[i] == 0) {
            continue;
        }

        ulong ntwid = ndx - two_partitions[i][0];
        ulong mtwid = ndy - two_partitions[i][1];
        ulong ptwid = ndz - two_partitions[i][2];

        total += multinomial_coef[i]*u_nmp_over_v(ntwid,mtwid,ptwid,x,y,z)*numerator_devs[i];
    }

    return total/(12*M_PI);
}

