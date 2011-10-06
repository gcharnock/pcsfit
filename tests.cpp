
#include "tests.hpp"
#include "model.hpp"
#include "model2.hpp"
#include "fit.hpp"

#include <boost/bind.hpp>
#include <iostream>

using namespace std;

RandomDist dist;

void check_derivative (PRNG prng,const Model* model) {

    //Chekcs the analytic and numerical derivates are equil
    for(unsigned long i = 0;i<10;i++) {
        double result;
        double* params             = (double*)alloca(model->size);
        double* gradient           = (double*)alloca(model->size);
		double* numerical_gradient = (double*)alloca(model->size);

        for(unsigned long j = 0; j< model->size; j++) {
            params[j] = dist(prng);
        }

		Vector3 evalAt(dist(prng),dist(prng),dist(prng));

        model->modelf(evalAt,params,&result,gradient);

        numerical_derivative(boost::bind(model->modelf_ND,evalAt,_1,_2),params,model->size,numerical_gradient);

        cout << "Result = " << endl;
        for(unsigned long i=0;i<model->size;i++) {
			cout << name_param(POINT_PARAM(i)) << " = " << params[i] <<
				" grad = " << gradient[i]
				 << " ceneral Differenceing gradient of " << numerical_gradient[i] << endl;
        }
        cout << "================================================================================" << endl;
    }
}

void check_error_derivate(PRNG prng,const Model* model,Multithreader<fdf_t>* pool) {
    //Chekcs the analytic and numerical derivates are equil
    for(unsigned long i = 0;i<10;i++) {
        double result;
        double* params             = (double*)alloca(model->size);
        double* gradient           = (double*)alloca(model->size);
		double* numerical_gradient = (double*)alloca(model->size);;

        for(unsigned long j = 0; j< model->size; j++) {
            params[j] = dist(prng);
        }

        Dataset dataset;

        random_data(prng,point_model,params,5*i,&dataset);

        ErrorContext context;
        context.dataset = &dataset;
        context.model   = model;
        context.params  = params;
        context.pool    = pool;

        eval_error(&context,&result,gradient);

        /*
          This section may be affected by a wierd compiler bug where
          the first argument of numerical_derivative ends up pointing
          to zero.
         */

        HasNDerivative curried_eval_error = boost::bind(eval_error,&context,_1,_2);
        numerical_derivative(&curried_eval_error,params,model->size,numerical_gradient);

        cout << "Result = " << endl;
        for(unsigned long i=0;i<model->size;i++) {
			cout << name_param(POINT_PARAM(i)) << " = " << params[i] <<
				" grad = " << gradient[i]
				 << " ceneral Differenceing gradient of " << numerical_gradient[i] << endl;
        }
        cout << "================================================================================" << endl;
    }
}

//Perform a sanity check on the models. Given teh same paramiter set
//the gaussian model should converged to the point model as stddev
//tends to 0
void testModel(long seed) {
	PRNG prng(seed);  //We should actually be recycling prng rather
					  //than the seed

	PointModel pm = PointModel::randomModel(seed+10);
    double pm2[8];

    pm2[PARAM_X]     = pm.metal.x;
    pm2[PARAM_Y]     = pm.metal.y;
    pm2[PARAM_Z]     = pm.metal.z;
           
    pm2[PARAM_CHI1]  = dist(prng);
    pm2[PARAM_CHI2]  = dist(prng);
    pm2[PARAM_CHIXY] = dist(prng);
    pm2[PARAM_CHIXZ] = dist(prng);
    pm2[PARAM_CHIYZ] = dist(prng);



	cout << "================================================================================" << endl;

	cout << "Checking the gaussian is normalised" << endl;
	cout << "TODO" << endl;

	cout << "================================================================================" << endl;
	cout << "Point Model = " << pm << endl;

	GaussModel gm1;

	gm1.ax = pm.ax;
	gm1.rh = pm.rh;

	gm1.metal = pm.metal;
	gm1.setEulerAngles(pm.angle_x,pm.angle_y,pm.angle_z);

	gm1.stddev = 1;

	GaussModel gm05 = gm1;
	gm05.stddev = 0.5;

	GaussModel gm01 = gm1;
	gm01.stddev = 0.1;

	GaussModel gm005 = gm1;
	gm005.stddev = 0.05;

	GaussModel gm001 = gm1;
	gm001.stddev = 0.01;

	GaussModel gm0005 = gm1;
	gm0005.stddev = 0.005;

	GaussModel gm0001 = gm1;
	gm0001.stddev = 0.001;

	GaussModel gm00005 = gm1;
	gm00005.stddev = 0.0005;


	cout << "================================================================================" << endl;
	for(unsigned long i = 0;i<0;i++) {
		double x = dist(prng);
		double y = dist(prng);
		double z = dist(prng);

		double dx = pm.metal.x - x;
		double dy = pm.metal.y - y;
		double dz = pm.metal.z - z;

		double r = sqrt(dx*dx + dy*dy + dz*dz);

		cout << "Point selected = ("<< x << "," << y << "," << z << ")" 
			 << "Metal-point distance is " << r << endl;

		cout << "Point Model = " << pm.eval(x,y,z) << endl;
		cout << "gm1 Model = " << gm1.eval(x,y,z) << endl;
		cout << "gm05 Model = " << gm05.eval(x,y,z) << endl;
		cout << "gm01 Model = " << gm01.eval(x,y,z) << endl;
		cout << "gm005 Model = " << gm005.eval(x,y,z) << endl;
		cout << "gm001 Model = " << gm001.eval(x,y,z) << endl;
		cout << "gm0005 Model = " << gm0005.eval(x,y,z) << endl;
		cout << "gm0001 Model = " << gm0001.eval(x,y,z) << endl;
		cout << "gm00005 Model = " << gm00005.eval(x,y,z) << endl;
		cout << "================================================================================" << endl;
	}

    cout << "Point Model 2 = "; for(unsigned long i=0;i<8;i++){cout << pm2[i] << " ";} cout << endl;

    cout << "================================================================================" << endl;

    Multithreader<fdf_t> pool;
    cout << "Evautating the error function for a perfect match (should be zero)" << endl;
    for(unsigned long i = 0;i<10;i++) {
        double* pm2 = (double*)alloca(8*sizeof(double));
        for(unsigned long j=0;j<8;j++) {
            pm2[j] = dist(prng);
        }
        Dataset dataset;
        random_data(prng,point_model,pm2,40,&dataset);
        
        ErrorContext context;
        context.dataset = &dataset;
        context.params  = pm2;
        context.model   = &point_model;
        context.pool = &pool;

        double error;
        double gradient[8];

        eval_error(&context,&error,gradient);

        cout << "run " << i << ", error = " << error << " (should be zero)" << endl;
        for(unsigned long j = 0;j < 8;j++) {cout << gradient[j] << " ";}
        cout << endl;
    }

    cout << "================================================================================" << endl;

    check_derivative (prng,&point_model);
    check_derivative (prng,&gaussian_model);

    cout << "Evaulating the analytic and numerical derivatives of the error functional" << endl;
    check_error_derivate(prng,&point_model,&pool);
    check_error_derivate(prng,&gaussian_model,&pool);

    cout << "================================================================================" << endl;
    for(unsigned long i = 1;i<0;i++) {
        double pm2[8];
        double model_to_fit[8];
        double model_to_fit2[8];
        double model_final[8];

        for(unsigned long j=0;j<8;j++) {
            pm2[j]          = dist(prng);
            model_to_fit2[j] = model_to_fit[j] = dist(prng)*100;
            model_final[j]  = 0;
        }

        Dataset dataset;

        cout << "Generating a distom molecule with " << 10*i << " spins" << endl;
        random_data(prng,point_model,pm2,5*i,&dataset);
        
        ErrorContext context;
        context.dataset = &dataset;
        context.params  = model_to_fit;
        context.model   = &point_model;
        context.pool    = &pool;

        double errorFinal = 0;

        do_fit_with_grad(&context,model_final,&errorFinal);
        for(unsigned long j = 0;j < 8; j++) {
            cout << name_param(POINT_PARAM(j)) << ": real = " << pm2[j] << " start = " << model_to_fit2[j]
                 << " final = " << model_final[j] << endl;
        }
        cout << "The final error was: " << errorFinal << endl;
        cout << "================================================================================" << endl;
    }
}


