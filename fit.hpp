
#include "data.hpp"
#include "model2.hpp"
#include "threads.hpp"

#include <utility>


typedef std::pair<double,std::vector<double> > fdf_t;

struct ErrorContext {
    const Dataset* dataset;

    const Model* model;
    double* params;

    Multithreader<fdf_t>* pool;
};

void do_fit_with_grad(ErrorContext* context,double* optModel);
void do_fit_with_grad(ErrorContext* context,double* optModel,double* finalError);

void eval_error(ErrorContext* context,double* value, double* gradient);

