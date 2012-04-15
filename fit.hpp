
#ifndef FIT_HPP
#define FIT_HPP

#include "data.hpp"
#include "model.hpp"
#include "threads.hpp"

#include <gsl/gsl_multimin.h>
#include <utility>


typedef std::pair<double,std::vector<double> > fdf_t;

struct ErrorContext {
    ErrorContext()
        : dataset(NULL),model(NULL),params(NULL),
          params_to_fix(NULL),pool(NULL),modelOptions(NULL) {
    }
    const Dataset* dataset;    //The dataset to fit too

    const Model* model;   //The model to use
    const double* params; //Starting paramiters

    const bool* params_to_fix;

    bool rescale;

    Multithreader<fdf_t>* pool;
    const ModelOptions* modelOptions;
};

typedef bool (*OnIterate) (const ErrorContext*,unsigned long,gsl_multimin_fdfminimizer*);

void do_fit_with_grad(const ErrorContext* context,
                      double* optModel,
                      double* finalError,
                      OnIterate onIterate,
                      Vals* final_vals = NULL);

void eval_error(const ErrorContext* context,
                const double* params,
                double* value,
                double* gradient,
                Vals* final_vals = NULL);

void numerical_error_derivative(const ErrorContext* context,double* value, double* gradient);


#endif
