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
#include "model2.hpp"
#include "data.hpp"
#include "vis.hpp"
#include "minimiser.hpp"
#include "tests.hpp"

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
            minimiser = new GradientMinimiser<M>(boost::bind(&NumericalExperiment::minf,this,_1),
												 boost::bind(&NumericalExperiment::withGrad,this,_1,_2,_3),
												 small_steps,starting_model);
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

        for(unsigned long i = 0;true;i++) {
            p = minimiser->iterate();
            M m = p.second;
            if(visualThread.fw) {
                visualThread.fw->setCalcVals(calcVals,m.metal,m.angle_x,m.angle_y,m.angle_z);
            }
            cout << p.first << "\t" << m << endl;
        }
        return p;
    }

	double withGrad(const M& thisModel,bool useGrad,vector<double>& g) {
		//Prepare the work
		std::vector<boost::function<double()> > work;

		std::vector<double> ax_grads; ax_grads.resize(nuclei.size());
		std::vector<double> rh_grads; rh_grads.resize(nuclei.size());

		for(unsigned long i = 0;i<nuclei.size();i++) {
			work.push_back(boost::bind(&M::grad,
									   &thisModel,nuclei[i].x,nuclei[i].y,nuclei[i].z,
									   &ax_grads[i],&rh_grads[i]));
		}
		Vals results;
		results.resize(work.size());

		//Execute in parallel
		pool->map(work,results);

		//Now reduce the results
		double total = 0;
		double ax_grad = 0;
		double rh_grad = 0;
	
		for(unsigned long i = 0;i<expVals.size();i++) {
			double diff = results[i] - expVals[i];
			total += diff*diff;

			ax_grad += 2*ax_grads[i]*diff;
			rh_grad += 2*rh_grads[i]*diff;
		}
		g[0] = ax_grad;
		g[1] = rh_grad;
		cout << "grad = " << ax_grad << " " << rh_grad << endl;

		std::vector<double> packed = M::pack(thisModel);
		std::vector<double> hVec = M::getSmallSteps();;

		calcVals = results;
        for(unsigned long i = 2;i < packed.size();i++) {
			std::vector<double> vprime = packed;
			
            double df_by_di = 0;
            double h = hVec[i];


			vprime[i] = packed[i]+h;
            df_by_di += minf(M::unpack(vprime));

            vprime[i] = packed[i]-h;
            df_by_di -= minf(M::unpack(vprime));

            df_by_di/=(2*h);

            g[i] = df_by_di;
        }
		cout << "Ngrad = " << g[0] << " " << g[1] << endl;


		//Record the results
		fout << total << "\t\t" << thisModel << endl;
		return total;
	}

	double minf(const M& thisModel) {
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
		fout << total << "\t\t" << thisModel << endl;
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





//When set to true, created threads should exit
bool threadExit = false;





int main(int argc,char** argv) {
	cout.width(20);
	cout.fill(' ');
	cout.precision(10);
	cout << fixed << showpos;

	//Parse the command line and decide what to do
	using namespace boost::program_options;
	variables_map variablesMap;
	string modelType;
    string method;
	string filename;
	try {
		options_description optDesc("Options");

		optDesc.add_options()
			("run-tests","")
			("help","print this help page and exit")
			("gui","Visualise the fitting process with VTK")
			("model",value<string>(&modelType)->default_value("point"),
			 "Select the model to use")
            ("method",value<string>(&method)->default_value("bfgs"),"")
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
		if(variablesMap.count("run-tests") > 0) {
			testModel(1000);
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

	double inital_ax = -20504.4;
	double inital_rh = -17337.4;

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
			p_exp(nuclei,expVals,pm,variablesMap.count("gui") > 0,variablesMap["method"].as<string>()=="bfgs");

		p_exp.minimise();
	} else if(variablesMap["model"].as<string>() == "gauss") {
		GaussModel gm; //The best point fit
		gm.ax = inital_ax;   
		gm.rh = inital_rh;
		gm.metal = inital_metal;
		gm.setEulerAngles(inital_angle_x, inital_angle_y, inital_angle_z);
		gm.stddev = 0.1;

		NumericalExperiment<GaussModel>
			g_exp(nuclei,expVals,gm,variablesMap.count("gui") > 0,variablesMap["method"].as<string>()=="bfgs");

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
