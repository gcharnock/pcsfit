
#include "model_common.hpp"

#include <cmath>

#include "integrate.hpp"
#include "pointdev.hpp"
#include "maths.hpp"

using namespace std;

//================================================================================//

int Integrand3(const double xx[],double ff[],int ncomp, void* void_userdata) {
    Userdata* userdata = (Userdata*)(void_userdata);

    const double* params = userdata->params;
    double stddev = abs(params[PARAM_STDDEV]);

    unsigned long point_size = userdata->point_model->size;

    //The dummy variable
    Vec3d xyz(xx);
	double r2 = xyz.r2();

    //The free variable
    Vec3d evalAt = userdata->evalAt;

    //The metal
    Vec3d metal = Vec3d(params[0],params[1],params[2]);

    //When eval - xyz = metal there is a singularity in the point
    //model. singularity is the vector that xyz would have to be to
    //trigger the singularity.
    Vec3d singularity =  evalAt - metal;
    double singularity_r2 = singularity.r2();

    //Potential optimisation: pull the calculation of these out of the
    //integral where they can be
    double a_coef = 1/(stddev*stddev);   assert(isfinite(a_coef));
    double normalizer = pow(M_PI*stddev*stddev,-1.5);

    //theExp0 and rho0 are the normalized and unnormalized value of
    //the electron density at the singularity respectivly.
    double theExp0 = exp(-a_coef*singularity_r2);
    double rho0 = normalizer*theExp0;

    //A vector pointing from the centre of 1/r^3 to the center of the gaussian
    Vec3d expand_around = xyz - singularity;

    /*
    double expand_xx = expand_around.x()*expand_around.x();
    double expand_yy = expand_around.y()*expand_around.y();
    double expand_zz = expand_around.z()*expand_around.z();

    double expand_xy = expand_around.x()*expand_around.y();
    double expand_xz = expand_around.x()*expand_around.z();
    double expand_yz = expand_around.y()*expand_around.z();
    */

    double expand_r2 = expand_around.r2();

    //Evauate the model
    double f;
    double* gradient = ncomp == 1 ? NULL: (double*)alloca(point_size*sizeof(double));

    userdata->point_model->modelf(userdata->evalAt - xyz,params,&f,gradient);

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
        }
    }

    assert(isfinite(rho));
    assert(isfinite(rho0));
    assert(isfinite(f));


    Vec3d rho0_grad = firstGaussDerivative(stddev,expand_around) * normalizer;

    /*double d2fdxx,d2fdyy,d2fdzz,d2fdxy,d2fdxz,d2fdyz;
    secondGaussDerivative(stddev,normalizer*theExp0,expand_around,
                          &d2fdxx,
                          &d2fdyy,
                          &d2fdzz,
                          &d2fdxy,
                          &d2fdxz,
                          &d2fdyz);*/

    double toSub         = rho0 + correction*expand_around.dot(rho0_grad);
    /*double spacial_toSub = correction*(expand_xx*d2fdxx +
                                       expand_yy*d2fdyy +
                                       expand_zz*d2fdzz +
                  
                                       2*expand_xy*d2fdxy +
                                       2*expand_xz*d2fdxz +
                                       2*expand_yz*d2fdyz) / 2;*/ // 2 factorial

    if(abs((rho-toSub)/rho) < 1e-10) {
        cout << " Warning, possible loss of precision, rho = "
             << rho << " rho0 " << rho0 << " toSub = " << toSub 
             << "  diff = " << (rho-toSub) << endl;
        //cout << " computed diff = " << normalizer*gaussian_error_term_one(r,a,stddev) << endl;
        for(int i = 0; i < ncomp-1; i++) {
            cout << "ff[" << i+1 << "] = " << (rho-toSub) * gradient[i] << endl;
        }
        cout << "expand_r2 = " << expand_r2 << endl;
    }
    ff[0] = (rho - toSub)*f;
    assert(isfinite(ff[0]));

    if(ncomp == 1) {
        //If we don't need a gradient, we can stop here.
        return 0;
    }



    for(unsigned long i=1;i<9;i++){
        ff[i] = (rho - toSub) * gradient[i-1];
        //if(i == 1) cout << sqrt(prime_to_singularity2) << ", " << ff[i] << endl;
    }

    ff[1] = 0;
    ff[2] = 0;
    ff[3] = 0;

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

        double sxxxx = eval_point_model_dev_xyz(params,4,0,0,evalAt);
        double syyyy = eval_point_model_dev_xyz(params,0,4,0,evalAt);
        double szzzz = eval_point_model_dev_xyz(params,0,0,4,evalAt);

        double sxxyy = eval_point_model_dev_xyz(params,2,2,0,evalAt);
        double syyzz = eval_point_model_dev_xyz(params,0,2,2,evalAt);
        double szzxx = eval_point_model_dev_xyz(params,2,0,2,evalAt);

        double s2 = stddev*stddev;
        double s5 = stddev*s2*s2;
        
        double sgradient = (6/24.0) * s5 * (
                                            3.0/4.0 * (sxxxx+syyyy+szzzz) +
                                            1.0/2.0 * (sxxyy+syyzz+szzxx)
                                            );
        double h = 0.0002;
        double val_p,val_m;
        double val_2p,val_2m;
        double val_3p,val_3m;
        double val_4p,val_4m;
        
        double mute_params[GAUSS_SIZE];
        userdata.params = mute_params;

        memcpy(mute_params,params,sizeof(double)*GAUSS_SIZE);

        mute_params[PARAM_STDDEV] = stddev + 4*h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_4p,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev + 3*h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_3p,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev + 2*h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_2p,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev + h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_p,(void*)&userdata);


        mute_params[PARAM_STDDEV] = stddev - h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_m,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev - 2*h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_2m,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev - 3*h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_3m,(void*)&userdata);

        mute_params[PARAM_STDDEV] = stddev - 4*h;
        cuhreIntegrate(Integrand3,&bounds,1,&val_4m,(void*)&userdata);

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
