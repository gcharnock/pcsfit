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
#include <boost/program_options.hpp>
#include <sstream>
#include <cassert>

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
    VisualThread(){
		fw = NULL;
	}
	~VisualThread() {
		if(fw) delete fw;
	}
	void start() {
		fw = new FittingWindow;
	}
    void operator()() {
		if(!fw) {
			cerr << "Fitting window not allocated" << endl;
			assert(false);
			return;
		}
        fw->mainLoop();
    }
    FittingWindow* fw;
private:
    VisualThread(const VisualThread&);
};

//When set to true, created threads should exit
bool threadExit = false;

double minf(GaussModel thisModel) {
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
	if(visualThread.fw) {
		visualThread.fw->setCalcVals(*(calcVals.rbegin()),m.metal,m.angle_x,m.angle_y,m.angle_z);
	}
}


int main(int argc,char** argv) {
	//Parse the command line and decide what to do
	using namespace boost::program_options;
	variables_map variablesMap;
	try {
		options_description optDesc("Options");

		optDesc.add_options()
			("help","print this help page and exit")
			("gui","Visualise the fitting process with VTK")
			("model","")
			("input-file","")
			("grid","perform a grid search");

		positional_options_description positional;
		positional.add("model",1);
		positional.add("input-file",2);

		try {
			store(command_line_parser(argc,argv).options(optDesc).positional(positional).run(),variablesMap);
		} catch(std::exception& e) {
			cerr << e.what() << endl;

			cout << "Usage:" << endl;
			cout << optDesc << endl;
			return 1;

		}
		if(variablesMap.count("help")) {
			cout << "Usage:" << endl;
			cout << optDesc << endl;
			return 0;
		}
		if(variablesMap.count("input-file") != 1) {
			cout << "Specify exactly one input file" << endl;
			return 0;
		}

	} catch(std::exception& e) {
		//I'm not sure if this is possible
		cerr << "An error occured when trying to parse the command line" << endl;
	}

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
	if(variablesMap.count("gui") > 0) {
		visualThread.start();
		visualThread.fw->setNuclei(nuclei);
		visualThread.fw->setExpVals(expVals);
		boost::thread boostVisualThread(boost::ref(visualThread));
	}

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
