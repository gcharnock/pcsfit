
#include "model_common.hpp"

#include <cmath>

#include "integrate.hpp"
#include "maths.hpp"

using namespace std;

//================================================================================//

int Integrand3(const double xx[],double ff[],int ncomp, void* void_userdata) {
    Userdata* userdata = (Userdata*)(void_userdata);

    const double* params = userdata->params;
    double stddev = abs(params[PARAM_STDDEV]);

    unsigned long point_size = userdata->point_model->size;

    double x = xx[0];
    double y = xx[1];
    double z = xx[2];

    //if(abs(z) < 0.00001) {
    //static long c = 0;
    //c++;
    //if(c % 1000 == 0) {
    //iout << x << "  " << y << "  " << z << endl;
    //}
    //}

	double r2 = x*x + y*y + z*z;


    //The location of the singularity (optimisation: take these
    //calculations out of the integral)
    double singularity_x = (userdata->evalAt.x() - params[0]);;
    double singularity_y = (userdata->evalAt.y() - params[1]);
    double singularity_z = (userdata->evalAt.z() - params[2]);

    double singularity_r2 =
        singularity_x*singularity_x +
        singularity_y*singularity_y +
        singularity_z*singularity_z;

    //Potential optimisation: pull the calculation of these out of the integral
    double a_coef = 1/(stddev*stddev);   assert(isfinite(a_coef));
    double normalizer = pow(M_PI*stddev*stddev,-1.5);
    double theExp0 = exp(-a_coef*(singularity_r2));
    double rho0 = normalizer*theExp0;

    double stddev2 = stddev*stddev;
    double stddev4 = stddev2*stddev2;
    double stddev6 = stddev4*stddev2;

    double sx = singularity_x;
    double sy = singularity_y;
    double sz = singularity_z;

    //First derivative
    double d_rho0_x = -2*sx/stddev2 * rho0;
    double d_rho0_y = -2*sy/stddev2 * rho0;
    double d_rho0_z = -2*sz/stddev2 * rho0;

    double d2_rho0_xx = (4*sx*sx - 2*stddev2)/stddev4 * rho0;
    double d2_rho0_yy = (4*sy*sy - 2*stddev2)/stddev4 * rho0;
    double d2_rho0_zz = (4*sz*sz - 2*stddev2)/stddev4 * rho0;

    //Second Derivative
    double d2_rho0_xy = (4*sx*sy)/stddev4 * rho0;
    double d2_rho0_xz = (4*sx*sz)/stddev4 * rho0;
    double d2_rho0_yz = (4*sy*sz)/stddev4 * rho0;

    //Third Derivative
    double d3_rho0_xxx = (12*sx*stddev2 - 8*sx*sx*sx)/stddev6 * rho0;
    double d3_rho0_yyy = (12*sy*stddev2 - 8*sy*sy*sy)/stddev6 * rho0;
    double d3_rho0_zzz = (12*sz*stddev2 - 8*sz*sz*sz)/stddev6 * rho0;
                      
    double d3_rho0_xxy = (4*y*(stddev2-2*sx*sx))/stddev6 * rho0;
    double d3_rho0_xxz = (4*z*(stddev2-2*sx*sx))/stddev6 * rho0;
                      
    double d3_rho0_yyx = (4*x*(stddev2-2*sy*sy))/stddev6 * rho0;
    double d3_rho0_yyz = (4*z*(stddev2-2*sy*sy))/stddev6 * rho0;
                      
    double d3_rho0_zzx = (4*x*(stddev2-2*sz*sz))/stddev6 * rho0;
    double d3_rho0_zzy = (4*y*(stddev2-2*sz*sz))/stddev6 * rho0;
                      
    double d3_rho0_xyz = -(8*sx*sy*sz)/stddev6 * rho0;

    //A vector pointing from the centre of 1/r^3 to the center of the gaussian
    double expand_around_x = x - singularity_x;
    double expand_around_y = y - singularity_y;
    double expand_around_z = z - singularity_z;
    
    double expand_r2 =
        expand_around_x*expand_around_x +
        expand_around_y*expand_around_y +
        expand_around_z*expand_around_z;

    //Evauate the model
    double f;
    double* gradient = ncomp == 1 ? NULL: (double*)alloca(point_size*sizeof(double));

    Vec3d x_minus_xprime(userdata->evalAt.x() - x,
                         userdata->evalAt.y() - y,
                         userdata->evalAt.z() - z);

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

    double toSub =  correction*(rho0 + expand_around_x*d_rho0_x
                                +      expand_around_y*d_rho0_y
                                +      expand_around_z*d_rho0_z
                                  
                                /*+      xprime_to_singularity*xprime_to_singularity*d2_rho0_xx
                                +      yprime_to_singularity*yprime_to_singularity*d2_rho0_yy
                                +      zprime_to_singularity*zprime_to_singularity*d2_rho0_zz

                                +      2*xprime_to_singularity*yprime_to_singularity*d2_rho0_xy
                                +      2*xprime_to_singularity*zprime_to_singularity*d2_rho0_xz
                                +      2*yprime_to_singularity*zprime_to_singularity*d2_rho0_yz*/);

    if( abs((rho-toSub)/rho) < 1e-10) {
        //cout << " Warning, possible loss of precision, rho = "
        //<< rho << " toSub = " << toSub << "  diff = " << (rho-toSub);
        Vec3d r(x,y,z);
        Vec3d a(singularity_x,singularity_y,singularity_z);

        //cout << " computed diff = " << normalizer*gaussian_error_term_one(r,a,stddev) << endl;
    } else {
        ff[0] = (rho-toSub)*f;
    }

    assert(isfinite(ff[0]));
    if(ncomp == 1) {
        //If we don't need a gradient, we can stop here.
        return 0;
    }

    for(unsigned long i=1;i<9;i++){
        ff[i] = (rho-toSub) * gradient[i-1];

        //if(i == 1) cout << sqrt(prime_to_singularity2) << ", " << ff[i] << endl;
    }
    for(int i=0;i<ncomp;i++){assert(isfinite(ff[i]));}
    return 0;
}


void eval_gaussian_series(Vec3d evalAt,const double* params,double* value, double* gradient) {
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

	cuhreIntegrate(Integrand3,&bounds,ncomp,integral,(void*)&userdata);

    *value = integral[0];

    if(gradient != NULL) {
        memcpy(gradient,integral+1,POINT_SIZE*sizeof(double));

        //We compute the gradient with respect to stddev via a series
        //expansion.
        /*
        double sxxxx = eval_point_model_dev_xyz(params,4,0,0,evalAt);
        double syyyy = eval_point_model_dev_xyz(params,0,4,0,evalAt);
        double szzzz = eval_point_model_dev_xyz(params,0,0,4,evalAt);

        double sxxyy = eval_point_model_dev_xyz(params,2,2,0,evalAt);
        double syyzz = eval_point_model_dev_xyz(params,0,2,2,evalAt);
        double szzxx = eval_point_model_dev_xyz(params,2,0,2,evalAt);

        double s2 = stddev*stddev;
        double s5 = stddev*s2*s2;
        
        gradient[PARAM_STDDEV] = (6/24.0) * s5 * (
                                                  3.0/4.0 * (sxxxx+syyyy+szzzz) +
                                                  1.0/2.0 * (sxxyy+syyzz+szzxx)
                                                  );
        *//*
        double h = 0.0002;
        double val_plus,val_minus;
        double val_2plus,val_2minus;
        
        double mute_params[GAUSS_SIZE];
        userdata.params = mute_params;

        memcpy(mute_params,params,sizeof(double)*GAUSS_SIZE);

        mute_params[PARAM_STDDEV] = stddev + 2*h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_2plus,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev + h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_plus,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev - h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_minus,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev - 2*h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_2minus,(void*)&userdata);
        

        gradient[PARAM_STDDEV] = (val_2minus - 8*val_minus + 8*val_plus - val_2plus)/(12*h);*/
        gradient[PARAM_STDDEV] = 0;
    }
}
