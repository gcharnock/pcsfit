
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

void eval_error(ErrorContext* context,double* value, double gradient[]);

