
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


fdf_t worker(const ErrorContext* context,const double* params,bool withGradient,unsigned long jobN) {
    assert(jobN < context->dataset->nuclei.size());
    unsigned long size = context->model->size;

    double value;
    double* gradient_buff = withGradient ? (double*)alloca(size*sizeof(double)) : NULL;
    vector<double> gradient;

    context->model->modelf(context->dataset->nuclei[jobN], params, &value, gradient_buff,context->modelOptions);

    if(withGradient) {
        for(unsigned long i = 0;i < size;i++) {
            gradient.push_back(gradient_buff[i]);
        }
    }

    return pair<double,vector<double> >(value,gradient);
}


//Evaulates the error functional. Effectivly does everything one
//iteration of the minimiser needs to do except
void eval_error(const ErrorContext* context,const double* params,double* value, double* gradient) {
    unsigned long size = context->model->size;
    //Prepare the work
    vector<boost::function<pair<double,vector<double> >()> > work;
    vector<fdf_t> results;

    for(unsigned long i = 0;i<context->dataset->nuclei.size();i++) {
        work.push_back(boost::bind(&worker,context,params,gradient,i));
    }

    //Map the work to the result
    context->pool->map(work,results);

    //Reduce the results
    *value = 0;
    if(gradient != NULL) {
        for(unsigned long i = 0; i<size;i++) {
            gradient[i] = 0.0;
        }
    }

    for(unsigned long i = 0; i<context->dataset->nuclei.size();i++) {
        double modelMinusexpVal = results[i].first - context->dataset->vals[i];
        double sq_error=modelMinusexpVal*modelMinusexpVal; 
        //Divide by the square of the shift so that large shifts don't
        //contribute excessily
        double sq_shift = 1;//(context->dataset->vals[i]*context->dataset->vals[i]);
        (*value) += sq_error/sq_shift; 
        //cout << "Spin " << i << " error = " << modelMinusexpVal*modelMinusexpVal << endl;
        if(gradient != NULL) {
            for(unsigned long j = 0; j < size; j++) {
                //The gradient of the error is twice the gradient of the
                //model - the expeimental value
                gradient[j] += 2 * results[i].second[j] * modelMinusexpVal / sq_shift;

                assert(isfinite(gradient[j]));
                assert(isfinite(gradient[j]*gradient[j]));
            }
        }
    }
}

void numerical_error_derivative(const ErrorContext* context,double* value, double* gradient) {
	unsigned long size = context->model->size;

	double* params_mutable = (double*)alloca(size * sizeof(double));
	memcpy(params_mutable,context->params,size*sizeof(double));

    for(unsigned long i = 0;i<context->model->size;i++) {
		double h     = abs(params_mutable[i]*0.00001);
		double result_plus,result_minus;

		params_mutable[i] = context->params[i] + h;
		eval_error(context,params_mutable,&result_plus, NULL);
		params_mutable[i] = context->params[i] - h;
		eval_error(context,params_mutable,&result_minus,NULL);

		params_mutable[i] = context->params[i];

		gradient[i] = (result_plus-result_minus)/(2*h);
        assert(isfinite(gradient[i]));
        assert(isfinite(gradient[i]*gradient[i]));
	}
    eval_error(context,context->params,value,NULL);
}



//Translates to and from GSL language
void eval_error_fdf(const gsl_vector* x, void* voidContext,double *f, gsl_vector *df) {
    const ErrorContext* context = (const ErrorContext*)(voidContext);

    unsigned long size = context->model->size;

    double* rescaled_x  = (double*)alloca(size *sizeof(double));
    double* gradient    = (double*)alloca(size *sizeof(double));

    //Rescale the gradient
    for(unsigned long i = 0; i<size; i++) {
        if(context->rescale) {
            rescaled_x[i] = gsl_vector_get(x,i) * context->params[i];
        } else {
            rescaled_x[i] = gsl_vector_get(x,i);
        }
        assert(isfinite(rescaled_x[i]));
    }

    //Test df for nullness and if it's null have the fdf function
    //write to throwaway memory.
    eval_error(context,rescaled_x,f, gradient);

    //Rescale the gradient (if it is needed)
    if(df != NULL) {
        for(unsigned long i = 0; i<size; i++) {
            double re_gradient = context->rescale ? gradient[i]*context->params[i] : gradient[i];
            assert(isfinite(re_gradient));
            assert(isfinite(re_gradient*re_gradient));
            gsl_vector_set(df,i,re_gradient);
        }
    }
}


double eval_error_f  (const gsl_vector* v, void* void_p) {
    double value;
    eval_error_fdf(v,void_p,&value,NULL);
    return value;
}

void   eval_error_df (const gsl_vector* v, void* void_p, gsl_vector *df) {
    double value;
    eval_error_fdf(v,void_p,&value,df);
    return;
}


void do_fit_with_grad(const ErrorContext* context,double* optModel,double* finalError,OnIterate onIterate) {
    assert(context != NULL);
    assert(optModel != NULL);


    unsigned long size = context->model->size;


    for(unsigned long i = 0; i < size; i++) {
        assert(isfinite(context->params[i]));
    }

    gsl_vector* gslModelVec = gsl_vector_alloc(size);
    if(context->rescale) { 
        //Use rescaling, the inital vector is all ones
        for(unsigned long i = 0; i<size; i++) {
            gsl_vector_set(gslModelVec,i,1.0);
        }
    } else {
        //We don't rescale, just copy
        memcpy(gsl_vector_ptr(gslModelVec,0), context->params ,size*sizeof(double));
    }

    //Setup the function to be minimised
    gsl_multimin_function_fdf minfunc;
    minfunc.f      = &eval_error_f;
    minfunc.df     = &eval_error_df;
    minfunc.fdf    = &eval_error_fdf;
    minfunc.n      = size;
    minfunc.params = (void*)(context);
    
    //Setup the minimiser
    gsl_multimin_fdfminimizer* gslmin;
    //gsl_multimin_fdfminimizer_steepest_descent
    //gsl_multimin_fdfminimizer_conjugate_fr
    //gsl_multimin_fdfminimizer_vector_bfgs2
    gslmin = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2,size);
    
    gsl_multimin_fdfminimizer_set(gslmin,&minfunc,gslModelVec,0.1,0.2);

    
	for(unsigned long i = 0;;i++) {
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

