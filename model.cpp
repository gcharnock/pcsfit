
#include "model.hpp"
#include "data.hpp"
#include "cuba.h"
#include <utility>
#include <cassert>
#include <cmath>

#include <iostream>

#define LN2

#define NDIM 3
#define NCOMP 1
#define EPSREL 1e-5
#define EPSABS 1e-5
#define VERBOSE 0
#define LAST 4
#define MINEVAL 2000
#define MAXEVAL 5000000

#define KEY 0

using namespace std;

typedef pair<const GaussModel*,double*> pairModelDouble;

PointModel::PointModel() {
}

PointModel::~PointModel() {
}


double PointModel::eval(double x,double y,double z) const {
	double gx = x - metal.x;
	double gy = y - metal.y;
	double gz = z - metal.z;

	double gxr = mat[0]*gx + mat[1]*gy + mat[2]*gz;
	double gyr = mat[3]*gx + mat[4]*gy + mat[5]*gz;
	double gzr = mat[6]*gx + mat[7]*gy + mat[8]*gz;

	double gx2 = gxr*gxr;
	double gy2 = gyr*gyr;
	double gz2 = gzr*gzr;

	double r2 = gx2 + gy2 + gz2;

	//Prevent division by zero errors
	if(r2 == 0) {
		return 0;
	}

	double r = sqrt(r2);
	double r5 = r2*r2*r;

	return (ax*(2*gz2 - gx2 - gy2) + rh*(3.0/2.0)*(gx2-gy2))/r5;
}


std::vector<double> packPointModel(const PointModel& m) {
	std::vector<double> vec;
    vec.push_back(m.ax     );
    vec.push_back(m.rh     );
                           
    vec.push_back(m.metal.x);
    vec.push_back(m.metal.y);
    vec.push_back(m.metal.z);
                           
    vec.push_back(m.angle_x);
    vec.push_back(m.angle_y);
    vec.push_back(m.angle_z);

    return vec;
}

PointModel unpackPointModel(const std::vector<double>& v) {
    PointModel m;
    m.ax      = v[0];
    m.rh      = v[1];
             
    m.metal.x = v[2];
    m.metal.y = v[3];
    m.metal.z = v[4];
             
	m.setEulerAngles(v[5],v[6],v[7]);

    return m;
}

//================================================================================//

GaussModel::GaussModel() {
	ax = 1;
	rh = 0; 

	metal.x = 0;
	metal.y = 0;
	metal.z = 0;

	setEulerAngles(0,0,0);

	stddev = 1;

	cube_x_min = -5; cube_x_max = 5;
	cube_y_min = -5; cube_y_max = 5;
	cube_z_min = -5; cube_z_max = 5;
}

GaussModel::~GaussModel() {

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

	//Prevent division by zero errors
	if(r2 == 0) {
		ff[0] = 0;
		return 0;
	}

	double r = sqrt(r2);
	double r5 = r2*r2*r;
    
    double stddev2 = this_->stddev*this_->stddev;
	double g = exp((-x*x-y*y-z*z)/(2*stddev2)) / sqrt(2*M_PI*stddev2);
	double f = (this_->ax*(2*gz2 - gx2 - gy2) + this_->rh*(3.0/2.0)*(gx2-gy2))/r5;


	double result = g * f;
	//Scale the result
	ff[0] = result * 
		(this_->cube_x_max - this_->cube_x_min)*
		(this_->cube_y_max - this_->cube_y_min)*
		(this_->cube_z_max - this_->cube_z_min);
    assert(isfinite(ff[0]));

	return 0;
}


double GaussModel::eval(double x,double y,double z) const {
	double nuclearLocation[3];
	nuclearLocation[0] = x;
	nuclearLocation[1] = y;
	nuclearLocation[2] = z;

	int nregions, neval, fail;
	double integral[1], error[1], prob[1];

	pairModelDouble p(this,nuclearLocation);
	Cuhre(NDIM, NCOMP, Integrand, (void*)&p,
		  EPSREL, EPSABS, VERBOSE | LAST,
		  MINEVAL, MAXEVAL, KEY,
		  &nregions, &neval, &fail, integral, error, prob);
    assert(!isinf(*integral));
	return *integral;
}

std::vector<double> packGaussModel(const GaussModel& m) {
	std::vector<double> vec;
    vec.push_back(m.ax     );
    vec.push_back(m.rh     );
                           
    vec.push_back(m.metal.x);
    vec.push_back(m.metal.y);
    vec.push_back(m.metal.z);
                           
    vec.push_back(m.angle_x);
    vec.push_back(m.angle_y);
    vec.push_back(m.angle_z);

    vec.push_back(m.stddev);
    return vec;
}

GaussModel unpackGaussModel(const std::vector<double>& v) {
    GaussModel m;
    m.ax      = v[0];
    m.rh      = v[1];
             
    m.metal.x = v[2];
    m.metal.y = v[3];
    m.metal.z = v[4];
             
	m.setEulerAngles(v[5],v[6],v[7]);

    m.stddev = v[8];
    return m;
}
