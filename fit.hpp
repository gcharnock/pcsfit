
#include "data.hpp"
#include "model2.hpp"
#include "threads.hpp"

#include <utility>


typedef std::pair<double,std::vector<double> > fdf_t;

struct ErrorContext {
    const Nuclei* nuclei;
    const Vals* expvals;
    
    double* model;
    ModelF modelf;
    unsigned long size;
    
    Multithreader<fdf_t>* pool;
};

void do_fit_with_grad(ErrorContext* context,double* optModel);
void do_fit_with_grad(ErrorContext* context,double* optModel,double* finalError);

void eval_error(ErrorContext* context,double* value, double gradient[]);

