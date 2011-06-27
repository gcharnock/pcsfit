#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <gsl/gsl_multimin.h>
#include <vector>
#include <boost/bind.hpp>

#include "foreach.hpp"
#include "model.hpp"
#include "data.hpp"
#include "threads.hpp"


using namespace std;

//Global State

GaussModel gaussModel;
Nuclei nuclei;
Vals expVals;
Vals calcVals;

double minf(const gsl_vector * v, void * mtVoid) {
	Multithreader<double>* mt = (Multithreader<double>*)mtVoid;
	//Unpack the vector
    gaussModel.ax = gsl_vector_get(v, 0);
    gaussModel.rh = gsl_vector_get(v, 1);
    gaussModel.metal.x = gsl_vector_get(v, 2);
    gaussModel.metal.y = gsl_vector_get(v, 3);
    gaussModel.metal.z = gsl_vector_get(v, 4);
    gaussModel.setEulerAngles(gsl_vector_get(v, 5),gsl_vector_get(v, 6),gsl_vector_get(v, 7));
    gaussModel.exponant = gsl_vector_get(v, 8);


    /*cout << "x = " << gaussModel.metal.x;
    cout << " y = " << gaussModel.metal.y;
    cout << " z = " << gaussModel.metal.z;
    cout << " ax = " << gaussModel.ax;
    cout << " rh = " << gaussModel.rh;
    cout << " ax = " << gaussModel.angle_x;
    cout << " ay = " << gaussModel.angle_y;
    cout << " az = " << gaussModel.angle_z;
    cout << " exp = " << gaussModel.exponant;
    cout << endl;*/

	std::vector<boost::function<double()> > work;
	work.resize(nuclei.size());


	//Prepare the work
	for(unsigned long i = 0;i<nuclei.size();i++) {
		work[i] = boost::bind(&GaussModel::eval,&gaussModel,nuclei[i].x,nuclei[i].y,nuclei[i].z,1e-4);
	}
	//Execute in parallel
	mt->map(work,calcVals);
	//Reduce
	double total = 0;
	foreach(double val,calcVals) {
		total += val;
	}
	return total;
 }


int main() {
	//Start the thread pool
	Multithreader<double> mt;
	//Setup the data

	//pair<Nuclei,Vals> pair_nv = fakeData(&gaussModel,100);
    pair<Nuclei,Vals> pair_nv = loadData("dataset_one.inp");
    nuclei = pair_nv.first;
    expVals = pair_nv.second;
	calcVals.resize(expVals.size());
	
	//Setup the inital model

	gaussModel.ax = 0;
	gaussModel.metal.x = (nuclei.xmax + nuclei.xmin)/2;
	gaussModel.metal.y = (nuclei.ymax + nuclei.ymin)/2;
	gaussModel.metal.z = (nuclei.zmax + nuclei.zmin)/2;

    double sizex = nuclei.xmax - nuclei.xmin;
    double sizey = nuclei.ymax - nuclei.ymin;
    double sizez = nuclei.zmax - nuclei.zmin;

	gaussModel.cube_x_min = -sizex*2; gaussModel.cube_x_max = sizex*2;
	gaussModel.cube_y_min = -sizey*2; gaussModel.cube_y_max = sizey*2;
	gaussModel.cube_z_min = -sizez*2; gaussModel.cube_z_max = sizez*2;


    gaussModel.setEulerAngles(0,0,0);

	//Initalise the minimizer
	const static unsigned long nParams = 9;

    gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,nParams);
    gsl_multimin_function minfunc;

    minfunc.f = &minf;
    minfunc.n = nParams;
    minfunc.params = &mt;

    gsl_vector *vec = gsl_vector_alloc (nParams);
    gsl_vector_set (vec, 0, gaussModel.ax); 
    gsl_vector_set (vec, 1, gaussModel.rh); 
    gsl_vector_set (vec, 2, gaussModel.metal.x);
    gsl_vector_set (vec, 3, gaussModel.metal.y);
    gsl_vector_set (vec, 4, gaussModel.metal.z);
    gsl_vector_set (vec, 5, gaussModel.angle_x);
    gsl_vector_set (vec, 6, gaussModel.angle_y);
    gsl_vector_set (vec, 7, gaussModel.angle_z);
    gsl_vector_set (vec, 8, gaussModel.exponant);


    gsl_vector *step_size = gsl_vector_alloc (nParams);
    gsl_vector_set (step_size, 0, 0.1);
    gsl_vector_set (step_size, 1, 0.1);
    gsl_vector_set (step_size, 2, 0.1);
    gsl_vector_set (step_size, 3, 0.1);
    gsl_vector_set (step_size, 4, 0.1);
    gsl_vector_set (step_size, 5, 0.1);
    gsl_vector_set (step_size, 6, 0.1);
    gsl_vector_set (step_size, 7, 0.1);
    gsl_vector_set (step_size, 8, 0.1);

    gsl_multimin_fminimizer_set (minimizer, &minfunc, vec, step_size);
	
	//Start the visualiser

	//Enter the main loop
	while(true) {
        gsl_multimin_fminimizer_iterate (minimizer);
	}
	return 0;
}
