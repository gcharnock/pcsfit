
#include "model.hpp"
#include "data.hpp"
#include "cuba.h"
#include <utility>
#include <math.h>

#include <iostream>

#define NDIM 3
#define NCOMP 1
#define EPSREL 1
#define EPSABS 1e-8
#define VERBOSE 0
#define LAST 4
#define MINEVAL 0
#define MAXEVAL 50000

#define KEY 0

using namespace std;

typedef pair<const GaussModel*,double*> pairModelDouble;

GaussModel::GaussModel() {
	ax = 1;
	rh = 0; 

	metal.x = 0;
	metal.y = 0;
	metal.z = 0;

	setEulerAngles(0,0,0);

	ax = 1;
	rh = 0;

	exponant = 1;

	cube_x_min = -5; cube_x_max = 5;
	cube_y_min = -5; cube_y_max = 5;
	cube_z_min = -5; cube_z_max = 5;
}

void GaussModel::setEulerAngles(double angle_x,double angle_y,double angle_z) {
	//From the matrix and quaternion FAQ
	double A       = cos(angle_x);
    double B       = sin(angle_x);
    double C       = cos(angle_y);
    double D       = sin(angle_y);
    double E       = cos(angle_z);
    double F       = sin(angle_z);
    double AD      = A * D;
    double BD      = B * D;
    mat[0]  =   C * E;
    mat[1]  =  -C * F;
    mat[2]  =   D;

    mat[3]  =  BD * E + A * F;
    mat[4]  = -BD * F + A * E;
    mat[5]  =  -B * C;
    mat[6]  = -AD * E + B * F;
    mat[7]  =  AD * F + B * E;
    mat[8] =   A * C;

	cout << mat[0] << mat[1] << mat[2] << endl;
	cout << mat[3] << mat[4] << mat[5] << endl;
	cout << mat[6] << mat[7] << mat[8] << endl;
}


//Since we need to pass this function to a C api it can't actually be
//a member function. We'll pass the this pointer explicitly to fake a
//member
int Integrand(const int *ndim, const double xx[],
			  const int *ncomp, double ff[], void *userdata) {
	pairModelDouble* p = (pairModelDouble*)userdata;
	const GaussModel* this_ = p->first;
	double* nuclearLocation = p->second;

	//Apply the transformation from the unit cube [0,1]^2. x,y,z is
	//the dummy variable in molecular coordinates
	double x = xx[0]*(this_->cube_x_max - this_->cube_x_min) + this_->cube_x_min;
	double y = xx[1]*(this_->cube_y_max - this_->cube_y_min) + this_->cube_y_min;
	double z = xx[2]*(this_->cube_z_max - this_->cube_z_min) + this_->cube_z_min;

	//The vector from the nucleous to the metal (we would use this as
	//the input to the point model). As we are convolving, subtract
	//the dummy variable
	double gx = (this_->metal.x - nuclearLocation[0]) - x;
	double gy = (this_->metal.y - nuclearLocation[1]) - y;
	double gz = (this_->metal.z - nuclearLocation[2]) - z;

	const double* mat = this_->mat;

	double gxr = mat[0]*gx + mat[1]*gy + mat[2]*gz;
	double gyr = mat[3]*gx + mat[4]*gy + mat[5]*gz;
	double gzr = mat[6]*gx + mat[7]*gy + mat[8]*gz;

	double gx2 = gxr*gxr;
	double gy2 = gyr*gyr;
	double gz2 = gzr*gzr;

	double r2 = gx2 + gy2 + gz2;

	if(r2 < this_->cutoff2) {
		ff[0] = 0;
		return 0;
	}

	double r = sqrt(r2);
	double r5 = r2*r2*r;

	double g = exp(-x*x-y*y-z*z);
	double f = (this_->ax*(2*gz2 - gx2 - gy2) + this_->rh*(3.0/2.0)*(gx2-gy2))/r5;


	double result = g * f;
	//Scale the result
	ff[0] = result * 
		(this_->cube_x_max - this_->cube_x_min)*
		(this_->cube_y_max - this_->cube_y_min)*
		(this_->cube_z_max - this_->cube_z_min);

	return 0;
}


double GaussModel::eval(double x,double y,double z,double epsAbs) const {
	double nuclearLocation[3];
	nuclearLocation[0] = x;
	nuclearLocation[1] = y;
	nuclearLocation[2] = z;

	int comp, nregions, neval, fail;
	double integral[1], error[1], prob[1];

	pairModelDouble p(this,nuclearLocation);
	Cuhre(NDIM, NCOMP, Integrand, (void*)&p,
		  EPSREL, epsAbs, VERBOSE | LAST,
		  MINEVAL, MAXEVAL, KEY,
		  &nregions, &neval, &fail, integral, error, prob);
	return *integral;
}

void GaussModel::bulkEval(const Nuclei& nuclei,Vals& vals) const {
	//Todo, replace with 
	for(long i=0;i<nuclei.size();i++) {
		vals[i] = eval(nuclei[i].x,nuclei[i].y,nuclei[i].z,1e-4);
		cout << vals[i] << endl;
	}
	vals.updateMinMax();
}


