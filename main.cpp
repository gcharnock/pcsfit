#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <gsl/gsl_multimin.h>
#include <vector>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <sstream>
#include <cassert>
#include <boost/function.hpp>
#include <limits>

#include "foreach.hpp"
#include "threads.hpp"
#include "model.hpp"
#include "data.hpp"
#include "vis.hpp"
#include "minimiser.hpp"

#define logMsg(x) flog << __FILE__ << "(" << __LINE__ << ")" << x << endl;

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
ofstream flog;

/*********************************************************************************
 * NumericalExperiment class
 *********************************************************************************/

//
template<typename M> //M stands for Model
class NumericalExperiment {
public:
	NumericalExperiment(const Nuclei& _nuclei,
						const Vals& _expVals,
                        const M& starting_model,
						bool withGUI,bool withGradient)
		:  nuclei(_nuclei), expVals(_expVals) {

        std::vector<double> small_steps = M::getSmallSteps();

        if(withGradient) {
            minimiser = new GradientMinimiser<M>(boost::bind(&NumericalExperiment::minf,this,_1),small_steps,starting_model);
        } else {
            minimiser = new SimplexMinimiser<M>(boost::bind(&NumericalExperiment::minf,this,_1),small_steps,starting_model);
        }

		calcVals.resize(expVals.size());

		//Start the visualisation thread
		if(withGUI) {
			visualThread.start();
			visualThread.fw->setNuclei(nuclei);
			visualThread.fw->setExpVals(expVals);
			boost::thread boostVisualThread(boost::ref(visualThread));
		}
	}

	std::pair<double,M> minimise() {
        std::pair<double,M> p;

        for(unsigned long i = 0;i<2000;i++) {
            p = minimiser->iterate();
            M m = p.second;
            if(visualThread.fw) {
                visualThread.fw->setCalcVals(calcVals,m.metal,m.angle_x,m.angle_y,m.angle_z);
            }
        }
        return p;
    }

	double minf(const M& thisModel) {
		cout << thisModel << endl;

		static int paused = 0;
		if(paused > 2) {
			//usleep(10000);
		}
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
			double diff = expVals[i] - results[i];
			total += diff*diff;
		}

		calcVals = results;

		//Record the results
		fout << total << "\t\t" <<  thisModel << endl;

		return total;
	}

    //Initalise the minimizer
	MinimiserBase<M>* minimiser;

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


//Perform a sanity check on the models. Given teh same paramiter set
//the gaussian model should converged to the point model as stddev
//tends to 0
void testModel(long seed) {
	PRNG prng(seed);  //We should actually be recycling prng rather
					  //than the seed
	RandomDist rand;

	PointModel pm = PointModel::randomModel(seed+10);

	cout << "================================================================================" << endl;

	cout << "Checking the gaussian is normalised" << endl;
	cout << "TODO" << endl;

	cout << "================================================================================" << endl;
	cout << "Point Model = " << pm << endl;

	GaussModel gm1;

	gm1.ax = pm.ax;
	gm1.rh = pm.rh;

	gm1.metal = pm.metal;
	gm1.setEulerAngles(pm.angle_x,pm.angle_y,pm.angle_z);

	gm1.stddev = 1;

	GaussModel gm05 = gm1;
	gm05.stddev = 0.5;

	GaussModel gm01 = gm1;
	gm01.stddev = 0.1;

	GaussModel gm005 = gm1;
	gm005.stddev = 0.05;

	GaussModel gm001 = gm1;
	gm001.stddev = 0.01;

	GaussModel gm0005 = gm1;
	gm0005.stddev = 0.005;

	GaussModel gm0001 = gm1;
	gm0001.stddev = 0.001;

	GaussModel gm00005 = gm1;
	gm00005.stddev = 0.0005;


	cout << "================================================================================" << endl;
	for(unsigned long i = 0;i<10;i++) {
		double x = rand(prng);
		double y = rand(prng);
		double z = rand(prng);

		double dx = pm.metal.x - x;
		double dy = pm.metal.y - y;
		double dz = pm.metal.z - z;

		double r = sqrt(dx*dx + dy*dy + dz*dz);

		cout << "Point selected = ("<< x << "," << y << "," << z << ")" 
			 << "Metal-point distance is " << r << endl;

		cout << "Point Model = " << pm.eval(x,y,z) << endl;
		cout << "gm1 Model = " << gm1.eval(x,y,z) << endl;
		cout << "gm05 Model = " << gm05.eval(x,y,z) << endl;
		cout << "gm01 Model = " << gm01.eval(x,y,z) << endl;
		cout << "gm005 Model = " << gm005.eval(x,y,z) << endl;
		cout << "gm001 Model = " << gm001.eval(x,y,z) << endl;
		cout << "gm0005 Model = " << gm0005.eval(x,y,z) << endl;
		cout << "gm0001 Model = " << gm0001.eval(x,y,z) << endl;
		cout << "gm00005 Model = " << gm00005.eval(x,y,z) << endl;
		cout << "================================================================================" << endl;
	}

}


//When set to true, created threads should exit
bool threadExit = false;





int main(int argc,char** argv) {
	//Parse the command line and decide what to do
	using namespace boost::program_options;
	variables_map variablesMap;
	string modelType;
    string method;
	string filename;
	try {
		options_description optDesc("Options");

		optDesc.add_options()
			("test-models","")
			("help","print this help page and exit")
			("gui","Visualise the fitting process with VTK")
			("model",value<string>(&modelType)->default_value("point"),
			 "Select the model to use")
            ("method",value<string>(&method)->default_value("gradient"),"")
			("random-model-type",
			 "Specify the type of the model to use for generating random data")
			("input-file",value<string>(&filename),"specify an input file")
			("random-data",
			 "Instead of loading a file, generate random data and try and fit that")
			("show-only","Show the data only, don't try and fit. Implies --gui");

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
		if(variablesMap.count("test-models") > 0) {
			testModel(7734);
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

	//Open the log file and params.log file
	flog.open("log.log");
	if(!flog.is_open()) {
		cerr << "Could not open log.log for writing" << endl;
		return 1;
	}

    //Start the thread pool
    pool = new Multithreader<double>;

	//Load the data
	pairNucVals data;
	PointModel randomPointModel;
	if(variablesMap.count("random-data") == 0) {
		logMsg("Opening the data");
		data = loadData(filename);
	} else {
		logMsg("Generating a random point model");
		randomPointModel = PointModel::randomModel(54321);
		logMsg("Model is:" << randomPointModel);
		data = fakeData(123456,randomPointModel,200);
    }
		
	Nuclei nuclei = data.first;
	Vals expVals = data.second;


	fout.open("params.log");
	if(!fout.is_open()) {
		cerr << "Could not open params.log for writing" << endl;
		return 1;
	}
	logMsg("Log file opened");

	//Set up the model

	double inital_ax = -3000;//-20504.4;
	double inital_rh = -1000;//-17337.4;

	Vector3 inital_metal = Vector3(4.165,18.875,17.180);
	double inital_angle_x = 6.46403;
	double inital_angle_y = 16.436;
	double inital_angle_z = 14.2135;

	if(modelType == "point") {
		PointModel pm;
		pm.ax = inital_ax;  
		pm.rh = inital_rh;
		pm.metal = inital_metal;
		pm.setEulerAngles(inital_angle_x, inital_angle_y, inital_angle_z);


		NumericalExperiment<PointModel>
			p_exp(nuclei,expVals,pm,variablesMap.count("gui") > 0,variablesMap.count("gradient") > 0);

		p_exp.minimise();
	} else if(variablesMap["model"].as<string>() == "gauss") {
		GaussModel gm; //The best point fit
		gm.ax = inital_ax;   
		gm.rh = inital_rh;
		gm.metal = inital_metal;
		gm.setEulerAngles(inital_angle_x, inital_angle_y, inital_angle_z);
		gm.stddev = 0.1;

		NumericalExperiment<GaussModel>
			g_exp(nuclei,expVals,gm,variablesMap.count("gui") > 0,variablesMap.count("gradient") > 0);

		g_exp.minimise();
	} else {
		cerr << "Unknown model typw" << endl;
		return 1;
	}

    
    //NB: This delete currently causes an assert error in boost, I
    //don't know why.
	delete pool;
    return 0;
}
