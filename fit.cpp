
#include "fit.hpp"
#include <cassert>
#include "threads.hpp"
#include <gsl/gsl_multimin.h>
#include <boost/function.hpp>

using namespace std;

/*
  Debuging BFGS:
  1) Check the scaling of the problem. Are all variables about 1.
  2) Check the size of the step the minimiser is taking
  3) Try other minimiser such as gradient decent and conjugulate gradients
 */

#define MAX_ITTER 80

fdf_t worker(const ErrorContext* context,const double* params,unsigned long jobN) {
    assert(jobN < context->dataset->nuclei.size());

    double value;
    vector<double> gradient;
    double* gradient_array = (double*)alloca(context->model->size * sizeof(double));

    context->model->modelf(context->dataset->nuclei[jobN], params, &value, gradient_array);
    
    for(unsigned long i = 0;i < context->model->size;i++){
        gradient.push_back(gradient_array[i]);
    }

    return fdf_t(value,gradient);
}


//Evaulates the error functional. Effectivly does everything one
//iteration of the minimiser needs to do except
void eval_error(const ErrorContext* context,const double* params,double* value, double* gradient) {
    //Prepare the work
    vector<boost::function<pair<double,vector<double> >()> > work;
    vector<fdf_t> results;

    for(unsigned long i = 0;i<context->dataset->nuclei.size();i++) {
        work.push_back(boost::bind(&worker,context,params,i));
    }

    //Map the work to the result
    context->pool->map(work,results);

    //Reduce the results
    *value = 0;
    for(unsigned long i = 0; i<context->model->size; i++) {
        gradient[i] = 0;
    }
    for(unsigned long i = 0; i<context->dataset->nuclei.size();i++) {
        double modelMinusexpVal = results[i].first - context->dataset->vals[i];
        (*value)+=modelMinusexpVal*modelMinusexpVal; //Squared error
        for(unsigned long j = 0; j < context->model->size; j++) {
            //The gradient of the error is twice the gradient of the
            //model - the expeimental value
            gradient[j] += 2 * results[i].second[j] * modelMinusexpVal;
        }
    }
}

void numerical_error_derivative(const ErrorContext* context,double* value, double* gradient) {
	unsigned long size = context->model->size;
    double* fake_gradient = (double*)alloca(size * sizeof(double));

	double* params_mutable = (double*)alloca(size * sizeof(double));
	memcpy(params_mutable,context->params,size*sizeof(double));

    for(unsigned long i = 0;i<context->model->size;i++) {
		double h     = abs(params_mutable[i]*0.00001);
		double result_plus,result_minus;

		params_mutable[i] = context->params[i] + h;
		eval_error(context,params_mutable,&result_plus,fake_gradient);
		params_mutable[i] = context->params[i] - h;
		eval_error(context,params_mutable,&result_minus,fake_gradient);

		params_mutable[i] = context->params[i];

		gradient[i] = (result_plus-result_minus)/(2*h);
	}
    eval_error(context,context->params,value,fake_gradient);
}



//Translates to and from GSL language
void eval_error_fdf(const gsl_vector* v, void* voidContext,double *f, gsl_vector *df) {
    double fake_gradient[MAX_PARAMS];

    ErrorContext* context = (ErrorContext*)(voidContext);

    //Test df for nullness and if it's null have the fdf function
    //write to throwaway memory.
    eval_error(context,gsl_vector_const_ptr(v,0),f, df==NULL ? fake_gradient : gsl_vector_ptr(df,0) );
}


double eval_error_f  (const gsl_vector* v, void* voidContext) {
    double value;
    eval_error_fdf(v,voidContext,&value,NULL);
    return value;
}

void   eval_error_df (const gsl_vector* v, void* voidContext, gsl_vector *df) {
    double value;
    eval_error_fdf(v,voidContext,&value,df);
    return;
}


void do_fit_with_grad(const ErrorContext* context,double* optModel,double* finalError,OnIterate onIterate) {
    assert(context != NULL);
    assert(optModel != NULL);

    unsigned long size = context->model->size;

    gsl_vector* gslModelVec = gsl_vector_alloc(size);
    memcpy(gsl_vector_ptr(gslModelVec,0), context->params ,size*sizeof(double));

    //Setup the function to be minimised
    gsl_multimin_function_fdf minfunc;
    minfunc.f      = &eval_error_f;
    minfunc.df     = &eval_error_df;
    minfunc.fdf    = &eval_error_fdf;
    minfunc.n      = size;
    minfunc.params = (void*)context;
    
    //Setup the minimiser
    gsl_multimin_fdfminimizer* gslmin;
    gslmin = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs,size);
    
    gsl_multimin_fdfminimizer_set(gslmin,&minfunc,gslModelVec,0.01,0.1);

    
    for(unsigned long i = 0; i < MAX_ITTER; i++) {
        if(gsl_multimin_fdfminimizer_iterate(gslmin) == GSL_ENOPROG) {
            break;
        }
		if(!onIterate(context,i,gslmin)) {
			break;
		}
    }

    gsl_vector* gslOptVector = gsl_multimin_fdfminimizer_x(gslmin);
    (*finalError) = gsl_multimin_fdfminimizer_minimum(gslmin);

    memcpy(optModel, gsl_vector_ptr(gslOptVector,0), size*sizeof(double));

    //Clean up
    gsl_multimin_fdfminimizer_free(gslmin);
	gsl_vector_free(gslModelVec);
}

