
#include "model_common.hpp"
#include "pointdev.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <cstring>
#include <cmath>
#include <cassert>
#include <alloca.h>
#include <fstream>
#include <limits>

using namespace std;


ofstream iout("int.dat",ios::out);



std::string name_param(int param) {
    switch(param) {
    case PARAM_X: return "x";
    case PARAM_Y: return "y";
    case PARAM_Z: return "z";
                 
    case PARAM_CHI1:  return "chi_1";
    case PARAM_CHI2:  return "chi_2";
    case PARAM_CHIXY: return "chi_xy";
    case PARAM_CHIXZ: return "chi_xz";
    case PARAM_CHIYZ: return "chi_yz";

    case PARAM_STDDEV: return "stddev";
    }
    assert(false);
    return "";
}

void numerical_derivative(Vec3d evalAt,const Model* model,const double* params,double* gradient) {
	double* params_mutable = (double*)alloca(model->size*sizeof(double));
	memcpy(params_mutable,params,model->size*sizeof(double));

    for(unsigned long i = 0;i<model->size;i++) {
		double h     = abs(params[i]*0.001);
		double result_plus,result_minus;


		params_mutable[i] = params[i] + h;
		model->modelf(evalAt,params_mutable,&result_plus ,NULL);

		params_mutable[i] = params[i] - h;
		model->modelf(evalAt,params_mutable,&result_minus,NULL);

		params_mutable[i] = params[i];

		gradient[i] = (result_plus-result_minus)/(2*h);
	}
}




int parse_params_file(const std::string& filename,const Model** model,std::vector<double>* params) {
    ifstream fin(filename.c_str());
    if(!fin.is_open()) {
        return PARAM_FILE_NOT_FOUND;
    }
    string model_name;
    fin >> model_name;
    if(model_name == "point") {
        *model = &point_model;
    } else if(model_name == "gauss_test") {
        cout << "Using gauss_test model " << endl;
        *model = &gaussian_model_testing;
    } else if(model_name == "gauss") {
        *model = &gaussian_model_series;
    } else {
        return  UNKNOWN_MODEL;
    }
    params->resize((*model)->size);
    for(unsigned long i = 0;i < (*model)->size; i++) {
        fin >> params->at(i);
        if(fin.eof()) {
            return NOT_ENOUGH_PARAMS;
        }
    }
    return PARSE_SUCESS;
}


void random_data(PRNG& prng,const Model& model,const double* params,unsigned long natoms,Dataset* dataset) {
	RandomDist rand;

    //We should center our imaginary spins around the dipole or
    //fitting will be very hard.
    double x = params[PARAM_X];
    double y = params[PARAM_Y];
    double z = params[PARAM_Z];
    
    dataset->nuclei.resize(natoms);
    dataset->vals.resize(natoms);

    for(unsigned long i = 0; i < natoms;i++) {
		dataset->nuclei[i] = Vec3d(x + 10*rand(prng),y + 10*rand(prng),z + 10*rand(prng));
        model.modelf(dataset->nuclei[i],params,&(dataset->vals[i]),NULL);
    }
}


//================================================================================//

void eval_gaussian_testing(Vec3d evalAt,const double* params,double* value, double* gradient) {
    if(gradient != NULL) {
        //Don't bother with an analytical gradient
        for(unsigned long i = 0; i < gaussian_model_testing.size;i++) {
            gradient[i] = 0;
        }
    }

    double normalizer = 0;
    double total      = 0;

    double stddev = params[PARAM_STDDEV];
    double a_coef = 1/(stddev*stddev);  
    /*
      Eval (2l+1)^3 times on a squard grid and spacings of stddev/4
     */

    long l = 5;

    for(long i = -l; i < l; i++) { // -5,-4, ... , 4,5 if l=5
        for(long j = -l; j < l; j++) {
            for(long k = -l; k < l; k++) {
                double f;

                double x = i * (stddev/2);
                double y = j * (stddev/2);
                double z = k * (stddev/2);
                
                double r = sqrt(x*x+y*y+z*z);
                
                if(r < stddev/4) {
                    continue;
                }

                double rho = exp(-a_coef*r*r);

                Vec3d x_minus_xprime(evalAt.x() - x,
                                     evalAt.y() - y,
                                     evalAt.z() - z);
                eval_point(x_minus_xprime,params,&f,NULL);

                total += f * rho;
                normalizer += rho;
            }
        }
    }
    *value = total / normalizer;
}


const Model gaussian_model = {eval_gaussian,9,"Gaussian Model"};
const Model gaussian_model_series = {eval_gaussian_series,9,"Gaussian Model Series"};
const Model gaussian_model_testing = {eval_gaussian_testing,9,"Gaussian Model"};
