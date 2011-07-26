
#include "model.hpp"
#include "cuba.h"
#include <utility>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>


#include <iostream>

#define LN2

#define NDIM 3
#define NCOMP 1
#define EPSREL 1e-5
#define EPSABS 1000
#define VERBOSE 0
#define LAST 4
#define MINEVAL 100000
#define MAXEVAL 5000000

#define KEY 0

using namespace std;

//Used to pass information to Intergral though the C code
struct IntegralDetails {
    const GaussModel* model;
    
    double x;
    double y;
    double z;
    
    double xmax;
    double xmin;
               
    double ymax;
    double ymin;
               
    double zmax;
    double zmin;
};

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

	return (   ax*(2*gz2 - gx2 - gy2) + rh*(3.0/2.0)*(gx2-gy2)  )/r5/M_PI/12;
}


std::vector<double> PointModel::pack(const PointModel& m) {
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

PointModel PointModel::unpack(const std::vector<double>& v) {
    PointModel m;
    m.ax      = v[0];
    m.rh      = v[1];
             
    m.metal.x = v[2];
    m.metal.y = v[3];
    m.metal.z = v[4];
             
	m.setEulerAngles(v[5],v[6],v[7]);

    return m;
}

ostream& operator <<(ostream& out,const PointModel& m) {
	out << m.ax       << " "
		<< m.rh       << " "
		<< m.metal.x  << " "
		<< m.metal.y  << " "
		<< m.metal.z  << " "
		<< m.angle_x  << " "
		<< m.angle_y  << " "
		<< m.angle_z;
	return out;
}

PointModel PointModel::randomModel(long seed) {
	PRNG prng(seed);
	RandomDist rand;
    
    PointModel m;
    m.ax = (1-2*rand(prng))*10000;
    m.rh = (1-2*rand(prng))*10000;

    m.metal.x = rand(prng);
    m.metal.y = rand(prng);
    m.metal.z = rand(prng);

    m.setEulerAngles(2*M_PI*rand(prng),M_PI*rand(prng),2*M_PI*rand(prng));

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
}

GaussModel::~GaussModel() {

}

//Since we need to pass this function to a C api it can't actually be
//a member function. We'll pass the this pointer explicitly to fake a
//member
int Integrand(const int *ndim, const double xx[],
			  const int *ncomp, double ff[], void *userdata) {
    IntegralDetails* intDet = (IntegralDetails*)userdata;
	const GaussModel* this_ = intDet->model;

	//Apply the transformation from the unit cube [0,1]^2. xp,yp,zp is
	//the dummy variable in molecular coordinates. The free variable
	//is intDet->xyz
	double xp = xx[0]*(intDet->xmax - intDet->xmin) + intDet->xmin;
	double yp = xx[1]*(intDet->ymax - intDet->ymin) + intDet->ymin;
	double zp = xx[2]*(intDet->zmax - intDet->zmin) + intDet->zmin;

	//Now we evaluate the point model at xp

	//Let gx be the vector from the metal to xp
	double gx = xp - this_->metal.x;
	double gy = yp - this_->metal.y;
	double gz = zp - this_->metal.z;

	const double* mat = this_->mat;

	//Rotate gx
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

	double f = (this_->ax*(2*gz2 - gx2 - gy2) + this_->rh*(3.0/2.0)*(gx2-gy2))/(12*M_PI*r5); assert(isfinite(f));

	//Now evaulate rho at x-xp where x is the free variable
	double xp_m_x = xp - intDet->x;
	double yp_m_y = yp - intDet->y;
	double zp_m_z = zp - intDet->z;

	assert(isfinite(xp_m_x));
	assert(isfinite(yp_m_y));
	assert(isfinite(zp_m_z));
	
    double a_coef = 1/(this_->stddev*this_->stddev);                           	assert(isfinite(a_coef));
	double g = exp(-a_coef*(xp_m_x*xp_m_x+
							yp_m_y*yp_m_y+
							zp_m_z*zp_m_z)) * pow(abs(a_coef)/M_PI,1.5);       	assert(isfinite(g));



	double result = g * f;
	//Scale the result
	ff[0] = result * 
		(intDet->xmax - intDet->xmin)*
		(intDet->ymax - intDet->ymin)*
		(intDet->zmax - intDet->zmin);
    assert(isfinite(ff[0]));

	return 0;
}


double GaussModel::eval(double x,double y,double z) const {
    IntegralDetails intDet;
    intDet.model = this;

    intDet.x = x;
    intDet.y = y;
    intDet.z = z;

    intDet.xmax = x + x+5*stddev;
	intDet.xmin = x - x-5*stddev;
    
    intDet.ymax = y + y+5*stddev;
    intDet.ymin = y - y-5*stddev;
    
    intDet.zmax = z + z+5*stddev;
    intDet.zmin = z - z-5*stddev;

	int nregions, neval, fail;
	double integral[1], error[1], prob[1];


	Cuhre(NDIM, NCOMP, Integrand, (void*)&intDet,
		  EPSREL, EPSABS, VERBOSE | LAST,
		  MINEVAL, MAXEVAL, KEY,
		  &nregions, &neval, &fail, integral, error, prob);

    assert(isfinite(*integral));

	return *integral;
}

std::vector<double> GaussModel::pack(const GaussModel& m) {
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

GaussModel GaussModel::unpack(const std::vector<double>& v) {
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
ostream& operator <<(ostream& out,const GaussModel& m) {
	out << m.ax       << " "
		<< m.rh       << " "
		<< m.metal.x  << " "
		<< m.metal.y  << " "
		<< m.metal.z  << " "
		<< m.angle_x  << " "
		<< m.angle_y  << " "
		<< m.angle_z  << " "
		<< m.stddev;
	return out;
}

GaussModel GaussModel::randomModel(long seed) {
	PRNG prng(seed);
	RandomDist rand;
    
    GaussModel m;
    m.ax = (1-2*rand(prng))*10000;
    m.rh = (1-2*rand(prng))*10000;

    m.metal.x = rand(prng);
    m.metal.y = rand(prng);
    m.metal.z = rand(prng);

    m.setEulerAngles(2*M_PI*rand(prng),M_PI*rand(prng),2*M_PI*rand(prng));

	m.stddev = 0.2*rand(prng);

    return m;
}

