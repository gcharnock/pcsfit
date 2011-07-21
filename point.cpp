

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
#include <limits>

#include "foreach.hpp"
#include "model.hpp"
#include "data.hpp"
#include "vis.hpp"
#include "minimiser.hpp"


using namespace std;
using namespace boost;

/*********************************************************************************
 * Global State
 *********************************************************************************/

//Nuclear Positions
Nuclei nuclei;

//Experimental Vals
Vals expVals;

//Calculated values to display
Vals calcVals;

PointModel lastModel;

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

double minf(PointModel thisModel) {

    //Now reduce the results
    double total = 0;
	
    Vals results;
    results.resize(nuclei.size());

    for(unsigned long i = 0;i<expVals.size();i++) {
        results[i] = thisModel.eval(nuclei[i].x,nuclei[i].y,nuclei[i].z);
        double diff = expVals[i] - results[i];
        total += diff*diff;
    }
    //Push the results into 
    calcVals = results;

	//Record the results
	fout << thisModel.ax       << " "
		 << thisModel.rh       << " "
		 << thisModel.metal.x  << " "
		 << thisModel.metal.y  << " "
		 << thisModel.metal.z  << " "
		 << thisModel.angle_x  << " "
		 << thisModel.angle_y  << " "
		 << thisModel.angle_z  << " " << endl;

    return total;
}

void onIterate(VisualThread& visualThread,const PointModel& m) {
    visualThread.fw.setCalcVals(calcVals,m.metal,m.angle_x,m.angle_y,m.angle_z);
}

int main() {

    //Load the data
	pairNucVals data = loadData("dataset_one.inp");
	
	nuclei = data.first;
	expVals = data.second;

	calcVals.resize(expVals.size());

	
	//Set up the model
	PointModel m;
	m.ax = 100.0;
	m.rh = 0.0;
	m.metal = Vector3((nuclei.xmin+nuclei.xmax)/2,
					   (nuclei.ymin+nuclei.ymax)/2,
					   (nuclei.zmin+nuclei.zmax)/2);
	m.setEulerAngles(0,0,0);

    //Start the visualisation thread
    VisualThread visualThread;
    visualThread.fw.setNuclei(nuclei);
    visualThread.fw.setExpVals(expVals);
    boost::thread boostVisualThread(boost::ref(visualThread));

    //Initalise the minimizer
	Minimiser<PointModel> minimiser(&unpackPointModel,
									&packPointModel,
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

    double best = numeric_limits<double>::infinity();
    PointModel bestModel;

    //Main loop
    for(long i = 0;i<10;i++) {
        for(long j = 0;j<10;j++) {
            for(long k = 0;k<10;k++) {
                double fx = i/10.0;
                double fy = j/10.0;
                double fz = k/10.0;
                m.metal.x = nuclei.xmin*(1-fx) + nuclei.xmax*fx;
                m.metal.y = nuclei.xmin*(1-fy) + nuclei.xmax*fy;
                m.metal.z = nuclei.xmin*(1-fz) + nuclei.xmax*fz;
                cout << "Starting: (" << m.metal.x << "," << m.metal.y << "," << m.metal.z << ")" << endl;
                
                std::pair<double,PointModel> thePair = minimiser.minimise(m);
                double found_min = thePair.first;
                if(best > found_min) {
                    bestModel = thePair.second;
                    best = found_min;

                    cout << "New Best Model:" << endl;
                    cout << "ax = " << bestModel.ax
                         << "rh = " << bestModel.rh

                         << "x = " << bestModel.metal.x
                         << "y = " << bestModel.metal.y
                         << "z = " << bestModel.metal.z

                         << "angle_x = " << bestModel.angle_x
                         << "angle_y = " << bestModel.angle_y
                         << "angle_z = " << bestModel.angle_z << endl << endl;

                }
                cout << "Finished, minimum = " << found_min << " best so far " << best <<endl;
            }
        }
    }
    return 0;
}
