
#include <string>
#include "data.hpp"

#define MAX_PARAMS 10

typedef void (*ModelF)    (Vector3 evalAt, double* model,double* value,double* gradient);
typedef void (*ModelF_ND) (Vector3 evalAt, double* model,double* value);

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


void numerical_derivative(Vector3 evalAt,double* model,ModelF modelf,unsigned long nparams,double * gradient);
