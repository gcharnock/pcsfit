
#ifndef MODEL2_HPP
#define MODEL2_HPP

#include <string>
#include "data.hpp"

#define MAX_PARAMS 10


typedef void (*ModelF)    (Vector3 evalAt, const double* model,double* value,double* gradient);

struct Model { //Encodes a function with an analytic gradient, and size paramiters
    ModelF    modelf;
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

enum AR_POINT_PARAM {
    PARAM_CHI_AX = PARAM_Z+1, 
    PARAM_CHI_RH, 
    PARAM_ALPHA,
    PARAM_BETA,
    PARAM_GAMMA
};


void axrh_to_cart(const double* axrh_params,double* cart_params);
void cart_to_axrh(const double* axrh_params,double* cart_params);

std::string name_param(POINT_PARAM param);

void eval_point(   Vector3 evalAt,const double* pm,double* value, double* gradient);
void eval_gaussian(   Vector3 evalAt,const double* model,double* value, double* gradient);

void random_data(PRNG& prng,const Model& model,const double* params,unsigned long natoms,Dataset* dataset);

//void numerical_derivative(Vector3 evalAt,const Model& model,unsigned long nparams,double * gradient);
void numerical_derivative(Vector3 evalAt,const Model* model,const double* params,double* gradient);

#endif

