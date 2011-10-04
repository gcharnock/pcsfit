
#include "fit.hpp"
#include <cassert>
#include "threads.hpp"
#include <gsl/gsl_multimin.h>
#include <boost/function.hpp>

using namespace std;



fdf_t worker(ErrorContext* context,unsigned long jobN) {
    assert(jobN < context->nuclei->size());

    double value;
    vector<double> gradient;
    double gradient_array[MAX_PARAMS];

    context->modelf(context->nuclei->at(jobN), context->model, &value, gradient_array);
    
    for(unsigned long i = 0;i < context->size;i++){
        gradient.push_back(gradient_array[i]);
    }

    return fdf_t(value,gradient);
}

void eval_error_fdf(const gsl_vector* v, void* voidContext,double *f, gsl_vector *df) {
    //ErrorContext* context = (ErrorContext*)voidContext;
}


double eval_error_f  (const gsl_vector* v, void* voidContext) {
    assert(false); //Hopefully we don't need this
    return 1.0;
}

void   eval_error_df (const gsl_vector* v, void* voidContext, gsl_vector *df) {
    assert(false); //Hopefully we don't need this
    return;
}

void eval_error(ErrorContext* context,double* value, double gradient[]) {
    //Prepare the work
    vector<boost::function<pair<double,vector<double> >()> > work;
    vector<fdf_t> results;

    for(unsigned long i = 0;i<context->nuclei->size();i++) {
        work.push_back(boost::bind(&worker,context,i));
    }

    //Map the work to the result
    context->pool->map(work,results);

    //Reduce the results
    *value = 0;
    for(unsigned long i = 0; i<context->size; i++) {
        gradient[i] = 0;
    }
    for(unsigned long i = 0; i<context->nuclei->size();i++) {
        (*value)+=(results[i].first - context->expvals->at(i));
        for(unsigned long j = 0; j < context->size; j++) {
            //The gradient of the error is twice the gradient of the
            //model times the error
            gradient[j] += 2 * results[i].second[j] * (*value);
        }
    }
}
/*
void do_fit() {
    gsl_vector* gslModelVec = gsl_vector_alloc(size);
    memcpy(gsl_vector_ptr(gslModelVec,0), model ,size*sizeof(double));

    //Setup the function to be minimised
    gsl_multimin_function_fdf minfunc;
    minfunc.f   = &eval_error_f;
    minfunc.df  = &eval_error_df;
    minfunc.fdf = &eval_error_fdf;
    minfunc.n   = size;
    minfunc.params = (void*)&context;
    
    //Setup the minimiser
    gsl_multimin_fminimizer* gslmin;
    gslmin = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,size);

    //Clean up
    gsl_multimin_fminimizer_free(gslmin);
}
*/
