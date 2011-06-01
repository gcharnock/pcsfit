
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

	ax = 1;
	rh = 0;

	exponant = 1;

	cube_x_min = -5; cube_x_max = 5;
	cube_y_min = -5; cube_y_max = 5;
	cube_z_min = -5; cube_z_max = 5;
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

	
	double gx2 = gx*gx;
	double gy2 = gy*gy;
	double gz2 = gz*gz;

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


