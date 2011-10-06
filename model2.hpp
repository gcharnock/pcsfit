
#ifndef MODEL2_HPP
#define MODEL2_HPP

#include <string>
#include "data.hpp"

#define MAX_PARAMS 10


typedef void (*ModelF)    (Vector3 evalAt, const double* model,double* value,double* gradient);
typedef void (*ModelF_ND) (Vector3 evalAt, const double* model,double* value);

struct Model { //Encodes a function with an analytic gradient, and size paramiters
    ModelF    modelf;
    ModelF_ND modelf_ND;
    unsigned long size;
    const char* name;
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

void eval_point(   Vector3 evalAt,const double* pm,double* value, double* gradient);
void eval_point_ND(Vector3 evalAt,const double* pm,double* value);

void eval_gaussian(   Vector3 evalAt,const double* model,double* value, double* gradient);
void eval_gaussian_ND(Vector3 evalAt,const double* model,double* value);

void random_data(PRNG& prng,const Model& model,const double* params,unsigned long natoms,Dataset* dataset);

//void numerical_derivative(Vector3 evalAt,const Model& model,unsigned long nparams,double * gradient);
void numerical_derivative(Vector3 evalAt,const Model* model,const double* params,double* gradient);

#endif

