
#ifndef FIT_HPP
#define FIT_HPP

#include "data.hpp"
#include "model2.hpp"
#include "threads.hpp"

#include <gsl/gsl_multimin.h>
#include <utility>


typedef std::pair<double,std::vector<double> > fdf_t;

struct ErrorContext {
    const Dataset* dataset; //The dataset to fit too

    const Model* model;   //The model to use
    const double* params; //Starting paramiters

    const bool* params_to_fix;

    bool rescale;

    Multithreader<fdf_t>* pool;
    const ModelOptions* modelOptions;
};

typedef bool (*OnIterate) (const ErrorContext*,unsigned long,gsl_multimin_fdfminimizer*);

void do_fit_with_grad(const ErrorContext* context,double* optModel,double* finalError,OnIterate onIterate);

void eval_error(const ErrorContext* context,const double* params,double* value, double* gradient);

void numerical_error_derivative(const ErrorContext* context,double* value, double* gradient);


#endif
