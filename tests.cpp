
#include "tests.hpp"
#include "model.hpp"
#include "model2.hpp"
#include "fit.hpp"

#include <boost/bind.hpp>
#include <iostream>

using namespace std;

RandomDist dist;

void check_minimum(PRNG& prng,const Model* model,Multithreader<fdf_t>* pool) {
    cout << "Evautating the error function for a perfect match (should be zero)" << endl;

    for(unsigned long i = 0;i<10;i++) {
        double* params = (double*)alloca(model->size*sizeof(double));
        for(unsigned long j=0;j<model->size;j++) {
            params[j] = dist(prng);
        }
        Dataset dataset;
        random_data(prng,*model,params,5,&dataset);
        
        ErrorContext context;
        context.dataset = &dataset;
        context.params  = params;
        context.model   = model;
        context.pool    = pool;
        context.rescale = false;

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

void check_derivative (PRNG& prng,const Model* model) {
    cout << "Checking the analytic and numerical derivaties match on " << model->name << endl;
    
    //Chekcs the analytic and numerical derivates are equil
    for(unsigned long i = 0;i<3;i++) {
        double result;
        double* params             = (double*)alloca(model->size*sizeof(double));
        double* gradient           = (double*)alloca(model->size*sizeof(double));
		double* numerical_gradient = (double*)alloca(model->size*sizeof(double));

        for(unsigned long j = 0; j< model->size; j++) {
            params[j] = dist(prng);
        }
        if(model == &gaussian_model) {
            params[PARAM_STDDEV] *= 0.01;
        }

		Vector3 evalAt(dist(prng),dist(prng),dist(prng));

        model->modelf(evalAt,params,&result,gradient);

        numerical_derivative(evalAt,model,params,numerical_gradient);

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
    cout << "Checking the analytic and numerical derivaties match for the error functional" << endl;

    //Chekcs the analytic and numerical derivates are equil
    for(unsigned long i = 1;i<4;i++) {
        double result;
        double* params             = (double*)alloca(model->size*sizeof(double));
        double* params2            = (double*)alloca(model->size*sizeof(double));
        double* gradient           = (double*)alloca(model->size*sizeof(double));
		double* numerical_gradient = (double*)alloca(model->size*sizeof(double));;

		//We need to generate two models, or else the gradient of the
		//error function will always be zero (and the numerical
		//gradient will be tiny and based off numerical errors)
        for(unsigned long j = 0; j< model->size; j++) {
            params[j] = dist(prng);
            params2[j] = dist(prng);
        }

        Dataset dataset;

		cout << "Creating a random molecule of " << 5*i << " spins" << endl;
        random_data(prng,point_model,params,5*i,&dataset);

        ErrorContext context;
        context.dataset = &dataset;
        context.model   = model;
        context.params  = params2;
        context.pool    = pool;

        eval_error(&context,params2,&result,gradient);

        numerical_error_derivative(&context,&result,numerical_gradient);

        cout << "Result = " << result << endl;
        for(unsigned long i=0;i<model->size;i++) {
			cout << name_param(POINT_PARAM(i)) << " = " << params[i];

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

void do_convergence(PRNG& prng,const Model* model,Multithreader<fdf_t>* pool) {
	unsigned long size = model->size;

    for(unsigned long i = 1;i<4;i++) {
        double* params_real        = (double*)alloca(size*sizeof(double));
        double* params_start       = (double*)alloca(size*sizeof(double));
        double* params_opt         = (double*)alloca(size*sizeof(double));

        for(unsigned long j=0;j<size;j++) {
            params_real[j]  = dist(prng);
            params_start[j] = params_real[j] + 0.1*dist(prng);
        }

        Dataset dataset;

        cout << "Generating a distom molecule with " << 10*i << " spins" << endl;
        random_data(prng,point_model,params_real,5*i,&dataset);

        cout << "Real params:" << endl;
        for(unsigned long j = 0;j<model->size;j++) {cout << params_real[j] << " ";}
        cout << endl;

        ErrorContext context;
        context.dataset = &dataset;
        context.params  = params_start;
        context.model   = model;
        context.pool    = pool;

		double errorStart = 666;
        double errorFinal = 666;
		double* gradient  = (double*)alloca(size*sizeof(double));

		eval_error(&context,params_start,&errorStart,gradient);

        do_fit_with_grad(&context,params_opt,&errorFinal,do_convergence_on_iterate);
        for(unsigned long j = 0;j < 8; j++) {
            cout << name_param(POINT_PARAM(j)) << ": real = " << params_real[j] << " start = " << params_start[j]
                 << " final = " << params_opt[j] << endl;
        }
		cout.setf(ios::floatfield,ios::scientific);

        cout << "The inital error was " << errorStart << " and the final error was: " << errorFinal << endl;

        cout << "================================================================================" << endl;
    }
}

void test_gaussian(PRNG& prng) {
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
            Vector3 pos = Vector3(rand(prng),rand(prng),rand(prng));
            nuclei[j] = pos;
            eval_gaussian(     pos,params,&(valsG[j])   ,NULL);
            eval_gaussian_testing(pos,params,&(valsTest[j]),NULL);
        }


        for(unsigned long j = 0; j<natoms; j++) {
            double x = nuclei[j].x;
            double y = nuclei[j].y;
            double z = nuclei[j].z;

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
void testModel(PRNG& prng,Multithreader<fdf_t>* pool) {
    cout << "================================================================================" << endl;

    //Check the gaussian model;
    //test_gaussian(prng);

    check_derivative (prng,&point_model);
    check_derivative (prng,&gaussian_model);

    //cout << "Evaulating the analytic and numerical derivatives of the error functional" << endl;
    check_error_derivate(prng,&point_model   ,pool);
    //check_error_derivate(prng,&gaussian_model,pool);

    check_minimum(prng,&point_model   ,pool);
    //check_minimum(prng,&gaussian_model,pool);
    
	do_convergence(prng,&point_model   ,pool);
	//do_convergence(prng,&gaussian_model,pool);
}


