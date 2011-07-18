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
#include <sstream>

#include "foreach.hpp"
#include "threads.hpp"
#include "model.hpp"
#include "data.hpp"
#include "vis.hpp"
#include "minimiser.hpp"


using namespace std;
using namespace boost;

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
 * Logging
 ********************************************************************************/

ofstream fout;


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

double minf(GaussModel thisModel) {
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
								   &thisModel,nuclei[i].x,nuclei[i].y,nuclei[i].z
								   ));
	}
	Vals results;
    results.resize(work.size());

	//Execute in parallel
    pool->map(work,results);

    //Now reduce the results
    double total = 0;
	
    for(unsigned long i = 0;i<expVals.size();i++) {
        double g = exp(-nuclei[i].x*nuclei[i].x
                       -nuclei[i].y*nuclei[i].y
                       -nuclei[i].z*nuclei[i].z);
        total+=g;
        if(expVals[i] > 1) {
            double diff = expVals[i] - results[i];
            total += diff*diff;
        }
    }
    //Push the results into 
    calcVals.push_back(results);
    models.push_back(thisModel);

	//Record the results
	fout << thisModel.ax       << " "
		 << thisModel.rh       << " "
		 << thisModel.metal.x  << " "
		 << thisModel.metal.y  << " "
		 << thisModel.metal.z  << " "
		 << thisModel.angle_x  << " "
		 << thisModel.angle_y  << " "
		 << thisModel.angle_z  << " "
		 << thisModel.stddev << endl;

	stringstream fnameStream;
	fnameStream << "results/points" << models.size() << ".log";
	ofstream fpoints(fnameStream.str().c_str());
	for(unsigned long i = 0;i<results.size();i++) {
		fpoints << expVals[i] << " " << results[i] << " " << (expVals[i]-results[i]) << endl;
	}

    return total;
}

void onIterate(VisualThread& visualThread,const GaussModel& m) {
    visualThread.fw.setCalcVals(*(calcVals.rbegin()),m.metal,m.angle_x,m.angle_y,m.angle_z);
}

int main() {
    //Start the thread pool
    pool = new Multithreader<double>;

    //Load the data
	pairNucVals data = loadData("dataset_one.inp");
	
	nuclei = data.first;
	expVals = data.second;

	calcVals.resize(expVals.size());

	
	//Set up the model
	GaussModel gm;
	gm.ax = 100.0;
	gm.rh = 0.0;
	gm.metal = Vector3((nuclei.xmin+nuclei.xmax)/2,
					   (nuclei.ymin+nuclei.ymax)/2,
					   (nuclei.zmin+nuclei.zmax)/2);
	gm.setEulerAngles(0,0,0);
	gm.stddev = 1;

    //Start the visualisation thread
    VisualThread visualThread;
    visualThread.fw.setNuclei(nuclei);
    visualThread.fw.setExpVals(expVals);
    boost::thread boostVisualThread(boost::ref(visualThread));

    //Initalise the minimizer
	Minimiser<GaussModel> minimiser(&unpackGaussModel,
									&packGaussModel,
									&minf,
									bind(onIterate,boost::ref(visualThread),_1));

	//Open the log file
	fout.open("params.log");
	if(!fout.is_open()) {
		cerr << "Could not open params.log for writing" << endl;
		return 1;
	}
	//Clear the results directory
	system("rm results/*");

    //Main loop
	minimiser.minimise(gm);

    delete pool;
    return 0;

}
