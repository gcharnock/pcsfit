
#ifndef MODEL2_HPP
#define MODEL2_HPP

#include <string>
#include "data.hpp"

#define MAX_PARAMS 10


struct ModelOptions {
    double absError;
    double relError;
};

typedef void (*ModelF) (Vec3d evalAt, const double* model, double* value,double* gradient,const ModelOptions* modelOptions);

struct Model { //Encodes a function with an analytic gradient, and size paramiters
    ModelF    modelf;
    unsigned long size;
    unsigned long var_size;
    unsigned long fixed_size;
    const char* name;
};


//Point models
extern const Model point_model;

//Extended models
extern const Model gaussian_model;
extern const Model gaussian_model_series;
extern const Model gaussian_model_testing; //A much simplier model that gaussian_model for testing
extern const Model gaussian_model_num_dev;
extern const Model gaussian_model_fix_s;

#define PARAM_X 0
#define PARAM_Y 1
#define PARAM_Z 2
                
#define PARAM_CHI1 3
#define PARAM_CHI2 4
#define PARAM_CHIXY 5
#define PARAM_CHIXZ 6
#define PARAM_CHIYZ 7

#define PARAM_STDDEV 8

#define POINT_SIZE 8
#define GAUSS_SIZE 9

enum {
    PARSE_SUCESS,
    UNKNOWN_MODEL,
    NOT_ENOUGH_PARAMS,
    PARAM_FILE_NOT_FOUND
};

//Parses a model file and obtains the type of model to use and the
//starting values. The second argument is a pointer to a const pointer
//to a model.
int parse_params_file(const std::string& filename,const Model** model,std::vector<double>* params);

void axrh_to_cart(const double* axrh_params,double* cart_params);
void cart_to_axrh(const double* axrh_params,double* cart_params);

std::string name_param(int param);

//Point models
void eval_point(   Vec3d evalAt,const double* params,double* value, double* gradient,const ModelOptions* modelOptions);

void eval_gaussian(Vec3d evalAt,const double* params,double* value, double* gradient,const ModelOptions* modelOptions);
void eval_gaussian_testing(Vec3d evalAt,const double* params,double* value, double* gradient,const ModelOptions* modelOptions);


void random_data(PRNG& prng,
                 const Model& model,
                 const double* params,
                 unsigned long natoms,
                 Dataset* dataset,
                 const ModelOptions* modelOptions);

//void numerical_derivative(Vector3 evalAt,const Model& model,unsigned long nparams,double * gradient);
void numerical_derivative(Vec3d evalAt,
                          const Model* model,
                          const double* params,
                          double* gradient,
                          const ModelOptions* modelOptions);

#endif

