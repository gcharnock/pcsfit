#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>
#include <unistd.h>
#include <gsl/gsl_multimin.h>
#include <vector>
#include <boost/thread.hpp>
#include <boost/ref.hpp>

#include "threads.hpp"
#include "model.hpp"
#include "data.hpp"
#include "vis.hpp"


using namespace std;

/*********************************************************************************
 * Global State
 *********************************************************************************/

//The thread pool
Multithreader<double>* pool = NULL;

//The current model
vector<GaussModel> models;

//Nuclear Positions
Nuclei nuclei;

//Experimental Vals
Vals expVals;

//Calculated values
vector<Vals> calcVals;

//The dimensions of the cube we will be integrating over.
double cube_x_min, cube_x_max;
double cube_y_min, cube_y_max;
double cube_z_min, cube_z_max;

/********************************************************************************
 * Visuliser
 *********************************************************************************/

struct VisualThread {
    VisualThread(){}
    void operator()() {
        fw.mainLoop();
    }
    FittingWindow fw;
private:
    VisualThread(const VisualThread&);
};

//When set to true, created threads should exit
bool threadExit = false;

double minf(const gsl_vector * v, void *) {
    //Unpack the model
    GaussModel thisModel;

    //Unpack the model
    thisModel.ax       = gsl_vector_get(v, 0);
    thisModel.rh       = gsl_vector_get(v, 1);
    thisModel.metal.x  = gsl_vector_get(v, 2);
    thisModel.metal.y  = gsl_vector_get(v, 3);
    thisModel.metal.z  = gsl_vector_get(v, 4);
    thisModel.setEulerAngles(gsl_vector_get(v, 5),
                             gsl_vector_get(v, 6),
                             gsl_vector_get(v, 7));
    thisModel.exponant = gsl_vector_get(v, 8);

    //TODO, we can reduce the integration region safley by quite a bit.
    thisModel.cube_x_min = thisModel.metal.x-100;
    thisModel.cube_x_max = thisModel.metal.x+100;
    thisModel.cube_y_min = thisModel.metal.y-100;
    thisModel.cube_y_max = thisModel.metal.y+100;
    thisModel.cube_z_min = thisModel.metal.z-100;
    thisModel.cube_z_max = thisModel.metal.z+100;

	//Prepare the work
	std::vector<boost::function<double()> > work;

	for(unsigned long i = 0;i<nuclei.size();i++) {
		work.push_back(boost::bind(
								   &GaussModel::eval,
								   &thisModel,nuclei[i].x,nuclei[i].y,nuclei[i].z,1e-4
								   ));
	}
	Vals results;
    results.resize(work.size());

	//Execute in parallel
    pool->map(work,results);

    //Now reduce the results
    double total = 0;
	
    for(unsigned long i = 0;i<expVals.size();i++) {
        double diff = expVals[i] - results[i];
        total += diff*diff;
    }
    //Push the results into 
    cout << "Pushing back results, results.size()=" << results.size() << endl;
    calcVals.push_back(results);
    models.push_back(thisModel);

    return total;
}

int main() {
    //Start the thread pool
    pool = new Multithreader<double>;

    //Load the data
	pairNucVals data = loadData("dataset_one.inp");
	
	nuclei = data.first;
	expVals = data.second;

	calcVals.resize(expVals.size());

    //Initalise the minimizer
    static const unsigned long nParams = 9;
    gsl_multimin_fminimizer* minimizer =
		gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,nParams);
    gsl_multimin_function minfunc;

    minfunc.f = &minf;
    minfunc.n = nParams;
    minfunc.params = NULL;

    gsl_vector *vec = gsl_vector_alloc (nParams);
    gsl_vector_set (vec, 0, 0.0); //Ax
    gsl_vector_set (vec, 1, 0.0); //Rh
    gsl_vector_set (vec, 2, (nuclei.xmin+nuclei.xmax)/2); //x
    gsl_vector_set (vec, 3, (nuclei.ymin+nuclei.ymax)/2); //y
    gsl_vector_set (vec, 4, (nuclei.zmin+nuclei.zmax)/2); //z
    gsl_vector_set (vec, 5, 0.0); //angle_x
    gsl_vector_set (vec, 6, 0.0); //angle_y
    gsl_vector_set (vec, 7, 0.0); //angle_z
    gsl_vector_set (vec, 8, 1.0); //exponant

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

    //Start the visualisation thread
    VisualThread visualThread;
    visualThread.fw.setNuclei(nuclei);
    visualThread.fw.setExpVals(expVals);
    boost::thread boostVisualThread(boost::ref(visualThread));


    //Main loop
     
    while(true) {
        gsl_multimin_fminimizer_iterate (minimizer);
        visualThread.fw.setCalcVals(*(calcVals.rbegin()),models.rbegin()->metal);
    }
    delete pool;
    return 0;
}
