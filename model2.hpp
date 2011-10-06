
#ifndef MODEL2_HPP
#define MODEL2_HPP

#include <string>
#include <boost/function.hpp>
#include "data.hpp"

#define MAX_PARAMS 10


typedef void (*ModelF)    (Vector3 evalAt, double* model,double* value,double* gradient);
typedef void (*ModelF_ND) (Vector3 evalAt, double* model,double* value);

//paramiters(in), value(out)
typedef boost::function<void(double*,double*)> HasNDerivative;

struct Model { //Encodes a function with an analytic gradient, and size paramiters
    ModelF    modelf;
    ModelF_ND modelf_ND;
    unsigned long size;
};

extern const Model point_model;
extern const Model gaussian_model;

enum POINT_PARAM {
    PARAM_X,    
    PARAM_Y,    
    PARAM_Z,    
                
    PARAM_CHI1, 
    PARAM_CHI2, 
    PARAM_CHIXY,
    PARAM_CHIXZ,
    PARAM_CHIYZ,

    PARAM_STDDEV
};

std::string name_param(POINT_PARAM param);

void eval_point(Vector3 evalAt, double pm[8],double* value, double gradient[8]);
void eval_point_ND(Vector3 evalAt,double pm[8],double* value);

void eval_gaussian(Vector3 evalAt,double* model,double* value, double gradient[9]);
void eval_gaussian_ND(Vector3 evalAt,double* model,double* value);

void random_data(PRNG prng,const Model& model,double* params,unsigned long natoms,Dataset* dataset);

//void numerical_derivative(Vector3 evalAt,const Model& model,unsigned long nparams,double * gradient);
void numerical_derivative(HasNDerivative f,double* params,unsigned long nparams,double* gradient);

#endif

