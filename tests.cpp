

#include "tests.hpp"
#include "model.hpp"
#include "fit.hpp"
#include "maths.hpp"
#include "foreach.hpp"

#include <boost/bind.hpp>
#include <vector>
#include <iostream>

using namespace std;
using namespace Eigen;

RandomDist dist;

void check_minimum(PRNG& prng,const Model* model,Multithreader<fdf_t>* pool,const ModelOptions* modelOptions) {
    cout << "Evautating the error function for a perfect match (should be zero)" << endl;

    for(unsigned long i = 0;i<10;i++) {
        double* params = (double*)alloca(model->size*sizeof(double));
        for(unsigned long j=0;j<model->size;j++) {
            params[j] = dist(prng);
        }
        Dataset dataset;
        random_data(prng,*model,params,5,&dataset,modelOptions);
        
        ErrorContext context;
        context.dataset = &dataset;
        context.params  = params;
        context.model   = model;
        context.pool    = pool;
        context.rescale = false;
        bool* params_to_fix = (bool*)alloca(model->size*sizeof(double));
        for(unsigned long j = 0; j < model->size; j++) {
            params_to_fix[j] = false;
        }
        context.params_to_fix = params_to_fix;


        double  error,n_error;
        double* gradient = (double*)alloca(model->size*sizeof(double));
        double* n_gradient = (double*)alloca(model->size*sizeof(double));

        eval_error(&context,params,&error,gradient);
        numerical_error_derivative(&context,&n_error,n_gradient);

        cout << "run " << i << ", error = " << error << ", n_error = " << n_error << " (should be zero)" << endl;
        cout << "Analytic: "; for(unsigned long j = 0;j < model->size;j++) {cout << gradient[j] << " ";} cout << endl;
        cout << "Numeric: "; for(unsigned long j = 0;j < model->size;j++) {cout << n_gradient[j] << " ";}  cout << endl;
    }
}

void check_derivative (PRNG& prng,const Model* model,const ModelOptions* modelOptions) {
    cout << "Checking the analytic and numerical derivaties match on " << model->name << endl;
    
    //Chekcs the analytic and numerical derivates are equil
    for(unsigned long i = 0;i<3;i++) {
        double result,point_result;
        double* params             = (double*)alloca(model->size*sizeof(double));
        double* gradient           = (double*)alloca(model->size*sizeof(double));
		double* numerical_gradient = (double*)alloca(model->size*sizeof(double));

        for(unsigned long j = 0; j< model->size; j++) {
            params[j] = dist(prng);
        }

		Vec3d evalAt(dist(prng),dist(prng),dist(prng));

        model->modelf(evalAt,params,&result,gradient,modelOptions);
        point_model.modelf(evalAt,params,&point_result,NULL,modelOptions);

        numerical_derivative(evalAt,model,params,numerical_gradient,modelOptions);

        double x = evalAt.x() - params[PARAM_X];
        double y = evalAt.y() - params[PARAM_Y];
        double z = evalAt.z() - params[PARAM_Z];

        cout << "evalAt = (" << evalAt.x() << "," << evalAt.y() << "," << evalAt.z() << ") r = "
             << sqrt(x*x+y*y+z*z) << endl;
        cout << "Result = " << result << " (point result = " << point_result <<  ")" << endl;
        for(unsigned long i=0;i<model->size;i++) {
			cout << name_param(i) << " = " << params[i] <<
				" grad = " << gradient[i]
				 << " ceneral Differenceing gradient of " << numerical_gradient[i] << endl;
        }
        cout << "================================================================================" << endl;
    }
}

void check_error_derivate(PRNG prng,const Model* model,Multithreader<fdf_t>* pool,const ModelOptions* modelOptions) {
    cout << "Checking the analytic and numerical derivaties match for the error functional" << endl;

    double* params             = (double*)alloca(model->size*sizeof(double));
    double* params2            = (double*)alloca(model->size*sizeof(double));
    double* gradient           = (double*)alloca(model->size*sizeof(double));
    double* numerical_gradient = (double*)alloca(model->size*sizeof(double));

    //Chekcs the analytic and numerical derivates are equil
    for(unsigned long i = 1;i<4;i++) {
        double result;

		//We need to generate two models, or else the gradient of the
		//error function will always be zero (and the numerical
		//gradient will be tiny and based off numerical errors)
        for(unsigned long j = 0; j< model->size; j++) {
            params[j] = dist(prng);
            params2[j] = dist(prng);
        }

        Dataset dataset;

		cout << "Creating a random molecule of " << 5*i << " spins" << endl;
        random_data(prng,point_model,params,5*i,&dataset,modelOptions);

        ErrorContext context;
        context.dataset = &dataset;
        context.model   = model;
        context.params  = params2;
        context.pool    = pool;
        bool* params_to_fix = (bool*)alloca(model->size*sizeof(double));
        for(unsigned long j = 0; j < model->size; j++) {
            params_to_fix[j] = false;
        }
        context.params_to_fix = params_to_fix;


        eval_error(&context,params2,&result,gradient);

        numerical_error_derivative(&context,&result,numerical_gradient);

        cout << "Result = " << result << endl;
        for(unsigned long i=0;i<model->size;i++) {
			cout << name_param(i) << " = " << params[i];

			cout.setf(ios::floatfield,ios::scientific);
			cout <<	" grad = " << gradient[i]
				 << " ceneral Differenceing gradient of " << numerical_gradient[i] << endl;
        }
        cout << "================================================================================" << endl;
    }
}

bool do_convergence_on_iterate(const ErrorContext* context,unsigned long i,gsl_multimin_fdfminimizer* minimizer) {
	double fx = gsl_multimin_fdfminimizer_minimum(minimizer);
	gsl_vector* g = gsl_multimin_fdfminimizer_gradient(minimizer);
    gsl_vector* x = gsl_multimin_fdfminimizer_x(minimizer);

	double norm = 0;
	for(unsigned long j = 0; j < g->size; j++) {
		norm += gsl_vector_get(g,j)*gsl_vector_get(g,j);
	}

    cout << i << "\t";
    for(unsigned long j = 0;j< x->size;j++) {cout << gsl_vector_get(x,j) << "\t";}
	cout << "\tf(x) = " << fx << "\t|grad| = " << norm << endl;

	return true;
}

void check_convergence(PRNG& prng,const Model* model,Multithreader<fdf_t>* pool,const ModelOptions* modelOptions) {
	unsigned long size = model->size;

    for(unsigned long i = 1;i<4;i++) {
        double* params_real        = (double*)alloca(size*sizeof(double));
        double* params_start       = (double*)alloca(size*sizeof(double));
        double* params_opt         = (double*)alloca(size*sizeof(double));

        for(unsigned long j=0;j<size;j++) {
            params_real[j]  = dist(prng);
            params_start[j] = params_real[j] + 0.5*dist(prng);
        }

        Dataset dataset;

        cout << "Generating a distom molecule with " << 10*i << " spins" << endl;
        random_data(prng,point_model,params_real,5*i,&dataset,modelOptions);

        cout << "Real params:" << endl;
        for(unsigned long j = 0;j<model->size;j++) {cout << params_real[j] << " ";}
        cout << endl;

        ErrorContext context;
        context.dataset       = &dataset;
        context.params        = params_start;
        context.model         = model;
        context.pool          = pool;
        context.rescale       = true;
        bool* params_to_fix = (bool*)alloca(size*sizeof(double));
        for(unsigned long j = 0; j < size; j++) {
            params_to_fix[j] = false;
        }
        context.params_to_fix = params_to_fix;
        

		double errorStart = 666;
        double errorFinal = 666;
		double* gradient  = (double*)alloca(size*sizeof(double));

		eval_error(&context,params_start,&errorStart,gradient);

        do_fit_with_grad(&context,params_opt,&errorFinal,do_convergence_on_iterate);
        for(unsigned long j = 0;j < 8; j++) {
            cout << name_param(j) << ": real = " << params_real[j] << " start = " << params_start[j];
            if(!context.rescale) {
                cout << " final = " << params_opt[j] << endl;
            } else {
                cout << " final = " << params_opt[j]*params_start[j] << endl;
            }
        }
		cout.setf(ios::floatfield,ios::scientific);

        cout << "The inital error was " << errorStart << " and the final error was: " << errorFinal << endl;

        cout << "================================================================================" << endl;
    }
}

void test_gaussian(PRNG& prng,const ModelOptions* modelOptions) {
	RandomDist rand;

    for(unsigned long i = 1; i<4; i++) {
        //Generate a random dataset using the Gaussian and Gaussian
        //test models. We expect them to match

        double* params = (double*)alloca(gaussian_model.size*sizeof(double));

        for(unsigned long j=0;j<gaussian_model.size;j++) {
            params[j]  = dist(prng);
        }
        unsigned long natoms = i*10;

        Nuclei nuclei;
        Vals valsG;
        Vals valsTest;

        nuclei.resize(natoms);
        valsG.resize(natoms);
        valsTest.resize(natoms);

        for(unsigned long j = 0; j < natoms;j++) {
            Vec3d pos = Vec3d(rand(prng),rand(prng),rand(prng));
            nuclei[j] = pos;
            eval_gaussian(        pos,params,&(valsG[j])   ,NULL,modelOptions);
            eval_gaussian_testing(pos,params,&(valsTest[j]),NULL,modelOptions);
        }


        for(unsigned long j = 0; j<natoms; j++) {
            double x = nuclei[j].x();
            double y = nuclei[j].y();
            double z = nuclei[j].z();

            double stddev = params[PARAM_STDDEV];
            double r_singularity = sqrt(x*x+y*y+z*z);

            cout << "(" << (r_singularity/stddev) << "," << valsG[j] << "," << valsTest[j] << ")" << endl;
        }
        cout << endl;
        cout << "================================================================================" << endl;
    }
}

//Perform a sanity check on the models. Given teh same paramiter set
//the gaussian model should converged to the point model as stddev
//tends to 0
void testModel(PRNG& prng,Multithreader<fdf_t>* pool,const ModelOptions* modelOptions) {
    cout << "================================================================================" << endl;

    //Check the gaussian model;
    //test_gaussian(prng);

    check_derivative (prng,&point_model,modelOptions);
    check_derivative (prng,&gaussian_model,modelOptions);

    //cout << "Evaulating the analytic and numerical derivatives of the error functional" << endl;
    check_error_derivate(prng,&point_model   ,pool,modelOptions);
    check_error_derivate(prng,&gaussian_model,pool,modelOptions);

    check_minimum(prng,&point_model   ,pool,modelOptions);
    check_minimum(prng,&gaussian_model,pool,modelOptions);
    
	check_convergence(prng,&point_model   ,pool,modelOptions);
	check_convergence(prng,&gaussian_model,pool,modelOptions);
}



double testPoly6(Vec3d r) {
    //6th order term of the taylor expansion of exp(-x2-y2-z2)
    double x = r[0];
    double y = r[1];
    double z = r[2];

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;

    double x4 = x2*x2;
    double y4 = y2*y2;
    double z4 = z2*z2;

    double x6 = x4*x2;
    double y6 = y4*y2;
    double z6 = z4*z2;
    return -1/6.0*x6 - 1/6.0*y6 - 1/2.0*y4*z2 - 1/2.0*y2*z4 - 1/6.0*z6
        - 1/2.0*(y2 +z2)*x4 - 1/2.0*(y4 + 2*y2*z2 + z4)*x2;
}

