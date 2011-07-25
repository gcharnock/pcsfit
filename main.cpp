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
#include <functional>

#include "foreach.hpp"
#include "threads.hpp"
#include "model.hpp"
#include "data.hpp"
#include "vis.hpp"
#include "minimiser.hpp"


using namespace std;
using namespace boost;

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

/*********************************************************************************
 * Multithreading
 *********************************************************************************/


//The thread pool
Multithreader<double>* pool = NULL;


/********************************************************************************
 * Logging
 ********************************************************************************/

ofstream fout;


/*********************************************************************************
 * NumericalExperiment class
 *********************************************************************************/

//


template<typename M> //M stands for Model
class NumericalExperiment {
public:
	NumericalExperiment(const Nuclei& _nuclei,
						const Vals& _expVals,
						bool withGUI,
						boost::function<M(const vector<double>&)> unpack,
						boost::function<vector<double>(const M&)> pack)
		: minimiser(unpack,
					pack,
					[=](const M& m){return this->minf(m);},
					[=](const M& m){this->onIterate(m);}),
		  nuclei(_nuclei),
		  expVals(_expVals) {


		calcVals.resize(expVals.size());

		//Start the visualisation thread
		if(withGUI) {
			visualThread.start();
			visualThread.fw->setNuclei(nuclei);
			visualThread.fw->setExpVals(expVals);
			boost::thread boostVisualThread(boost::ref(visualThread));
		}
	}

	std::pair<double,M> minimise(const M& input) {
		return  minimiser.minimise(input);
	}
	double minf(const M& thisModel) {
		static int paused = 0;
		if(paused > 2) sleep(1);
		paused++;

		//Prepare the work
		std::vector<boost::function<double()> > work;

		for(unsigned long i = 0;i<nuclei.size();i++) {
			work.push_back(boost::bind(
									   &M::eval,
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
			if(expVals[i] > 1) {
				double diff = expVals[i] - results[i];
				total += diff*diff;
			}
		}

		calcVals = results;

		//Record the results
		fout << total << "\t\t" <<  thisModel << endl;

		return total;
	}
	void onIterate(const M& m) {
		if(visualThread.fw) {
			visualThread.fw->setCalcVals(calcVals,m.metal,m.angle_x,m.angle_y,m.angle_z);
		}
	}

    //Initalise the minimizer
	Minimiser<M> minimiser;

	//Nuclear Positions
	Nuclei nuclei;

	//Experimental Vals
	Vals expVals;

	//Calculated values
	Vals calcVals;


	//The visualiser
	VisualThread visualThread;
private:
	NumericalExperiment(const NumericalExperiment&) {};
};



//When set to true, created threads should exit
bool threadExit = false;



int main(int argc,char** argv) {
	//Parse the command line and decide what to do
	using namespace boost::program_options;
	variables_map variablesMap;
	string modelType;
	string filename;
	try {
		options_description optDesc("Options");

		optDesc.add_options()
			("help","print this help page and exit")
			("gui","Visualise the fitting process with VTK")
			("model",value<string>(&modelType)->default_value("point"),"Select the model to use")
			("input-file",value<string>(&filename),"specify an input file")
			("random-data","Instead of loading a file, generate random data and try and fit that")
			("grid","perform a grid search");

		try {
			store(command_line_parser(argc,argv).options(optDesc).run(),variablesMap);
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
		if(variablesMap.count("input-file") != 1 && 
		   variablesMap.count("random-data") == 0) {
			cout << "Specify exactly one input file" << endl;
			return 0;
		}

	} catch(std::exception& e) {
		//I'm not sure if this is possible
		cerr << "An error occured when trying to parse the command line" << endl;
	}
	notify(variablesMap);

    //Start the thread pool
    pool = new Multithreader<double>;

	//Load the data
	pairNucVals data;
	if(variablesMap.count("random-data") == 0) {
		data = loadData(filename);
	} else {
		data = fakeData(123456,PointModel::randomModel(54321),200);
    }
		
	Nuclei nuclei = data.first;
	Vals expVals = data.second;

	//Open the log file
	fout.open("params.log");
	if(!fout.is_open()) {
		cerr << "Could not open params.log for writing" << endl;
		return 1;
	}

	//Set up the model

	if(modelType == "point") {
		PointModel pm;
		pm.ax = -5889.0;  
		pm.rh =  -5491.0;
		pm.metal = Vector3(4.165,18.875,17.180);/*Vector3((nuclei.xmin+nuclei.xmax)/2,
						   (nuclei.ymin+nuclei.ymax)/2,
						   (nuclei.zmin+nuclei.zmax)/2);*/
		pm.setEulerAngles(4.26073, 1.18864, -3.54324);

		std::function<PointModel(vector<double>)> 
			unpackPoint = [](vector<double> v) {
			PointModel m;
			m.ax      = v[0];
			m.rh      = v[1];
			
			m.metal.x = 4.165;
			m.metal.y = 18.875;
			m.metal.z = 17.180;
												  
			m.setEulerAngles(v[2],v[3],v[4]);
			
			return m;
		};

		std::function<vector<double>(PointModel)>
			packPoint = [](const PointModel& m){
			vector<double> vec;
			vec.push_back(m.ax     );
			vec.push_back(m.rh     );
			
			vec.push_back(m.angle_x);
			vec.push_back(m.angle_y);
			vec.push_back(m.angle_z);
			
			return vec;
		};

		NumericalExperiment<PointModel>
			p_exp(nuclei,expVals,variablesMap.count("gui") > 0,unpackPoint,packPoint);

		p_exp.minimise(pm);
	} else if(variablesMap["model"].as<string>() == "gauss") {
		GaussModel gm;
		gm.ax = -5889.0;
		gm.rh = -5491.0;
		gm.metal = Vector3(4.165,18.875,17.180);/*(nuclei.xmin+nuclei.xmax)/2,
						   (nuclei.ymin+nuclei.ymax)/2,
						   (nuclei.zmin+nuclei.zmax)/2);*/
		gm.setEulerAngles(4.26073, 1.18864, -3.54324);
		gm.stddev = 1;

		std::function<GaussModel(vector<double>)>
			packGauss = [](vector<double> v){
			GaussModel m;
			m.ax      = v[0];
			m.rh      = v[1];
			
			m.metal.x = 4.165;
			m.metal.y = 18.875;
			m.metal.z = 17.180;
												  
			m.setEulerAngles(v[2],v[3],v[4]);
												  
			m.stddev = v[5];
			return m;
		};

		std::function<vector<double>(const GaussModel&)>
			unpackGauss = [](const GaussModel& m){
			std::vector<double> vec;
			vec.push_back(m.ax     );
			vec.push_back(m.rh     );
													
			vec.push_back(m.angle_x);
			vec.push_back(m.angle_y);
			vec.push_back(m.angle_z);
													
			vec.push_back(m.stddev);
			return vec;
		};

		NumericalExperiment<GaussModel>
			g_exp(nuclei,expVals,variablesMap.count("gui") > 0,packGauss,unpackGauss);

		g_exp.minimise(gm);
	} else {
		cerr << "Unknown model typw" << endl;
		return 1;
	}

    
    //NB: This delete currently causes an assert error in boost, I
    //don't know why.
	delete pool;
    return 0;
}
