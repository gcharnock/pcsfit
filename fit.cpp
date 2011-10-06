
#include "fit.hpp"
#include <cassert>
#include "threads.hpp"
#include <gsl/gsl_multimin.h>
#include <boost/function.hpp>

using namespace std;

#define MAX_ITTER 40

fdf_t worker(ErrorContext* context,unsigned long jobN) {
    assert(jobN < context->dataset->nuclei.size());

    double value;
    vector<double> gradient;
    double gradient_array[MAX_PARAMS];

    context->model->modelf(context->dataset->nuclei[jobN], context->params, &value, gradient_array);
    
    for(unsigned long i = 0;i < context->model->size;i++){
        gradient.push_back(gradient_array[i]);
    }

    return fdf_t(value,gradient);
}


//Evaulates the error functional. Effectivly does everything one
//iteration of the minimiser needs to do except
void eval_error(ErrorContext* context,double* value, double* gradient) {
    //Prepare the work
    vector<boost::function<pair<double,vector<double> >()> > work;
    vector<fdf_t> results;

    for(unsigned long i = 0;i<context->dataset->nuclei.size();i++) {
        work.push_back(boost::bind(&worker,context,i));
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
    cout.setf(ios::floatfield,ios::scientific);
    cout << *value << endl;
}

//Translates to and from GSL language
void eval_error_fdf(const gsl_vector* v, void* voidContext,double *f, gsl_vector *df) {
    double fake_gradient[MAX_PARAMS];

    ErrorContext* context = (ErrorContext*)(voidContext);
    memcpy(context->params,gsl_vector_const_ptr(v,0),context->model->size*sizeof(double));

    //Test df for nullness and if it's null have the fdf function
    //write to throwaway memory.
    eval_error(context,f, df==NULL ? fake_gradient : gsl_vector_ptr(df,0) );
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


void do_fit_with_grad(ErrorContext* context,double* optModel,double* finalError) {
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
    gslmin = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2,size);
    
    gsl_multimin_fdfminimizer_set(gslmin,&minfunc,gslModelVec,1,0.1);

    
    for(unsigned long i = 0; i < MAX_ITTER; i++) {
        if(gsl_multimin_fdfminimizer_iterate(gslmin) == GSL_ENOPROG) {
            break;
        }
        if(gsl_multimin_test_gradient(gsl_multimin_fdfminimizer_gradient(gslmin), 1e-20) == GSL_SUCCESS) {
            break;
        }
    }

    gsl_vector* gslOptVector = gsl_multimin_fdfminimizer_x(gslmin);
    (*finalError) = gsl_multimin_fdfminimizer_minimum(gslmin);

    memcpy(optModel, gsl_vector_ptr(gslOptVector,0), size*sizeof(double));

    //Clean up
    gsl_multimin_fdfminimizer_free(gslmin);
}

