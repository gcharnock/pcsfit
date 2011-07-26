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
#include <functional>
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
						bool withGUI,
						std::function<M(const vector<double>&)> unpack,
						std::function<vector<double>(const M&)> pack)
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

	std::pair<double,M> multiMinimise(const M& input,unsigned long seed,std::function<M(unsigned long)> makeRandom) {
		double best = numeric_limits<double>::infinity();
		M bestModel;
		for(unsigned long i = 0;i<1000;i++) {
			//fout << "------------------------------------------------------------" << endl;
			M thisModel = makeRandom(seed+i);


			auto output = minimiser.minimise(thisModel);

			cout << thisModel << endl;
			cout << output.second << endl;
			cout << endl;


			double error = output.first;
			if(error < best) {
				logMsg("Found a new better model, with error " << error);
				logMsg("Starting model was " << thisModel);
				logMsg("Final model was " << output.second);
				best = error;
				bestModel = output.second;
			}
			//fout << endl;
			//fout << endl;
			//fout << endl;
		}
		return pair<double,M>(best,bestModel);
	}

	double minf(const M& thisModel) {
		cout << thisModel << endl;

		static int paused = 0;
		if(paused > 2) {
			usleep(10000);
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


//Perform a sanity check on the models. Given teh same paramiter set
//the gaussian model should converged to the point model as stddev
//tends to 0
void testModel(long seed) {
	PRNG prng(seed);  //We should actually be recycling prng rather
					  //than the seed
	RandomDist rand;

	PointModel pm = PointModel::randomModel(seed+10);

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

		cout << "Point selected = ("<< x << "," << y << "," << z << ")" << endl;

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
	string filename;
	try {
		options_description optDesc("Options");

		optDesc.add_options()
			("test-models","")
			("help","print this help page and exit")
			("gui","Visualise the fitting process with VTK")
			("model",value<string>(&modelType)->default_value("point"),
			 "Select the model to use")
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

	if(modelType == "point") {
		PointModel pm;
		pm.ax = -5889.0;  
		pm.rh =  -5491.0;
		pm.metal = Vector3(4.165,18.875,17.180);/*Vector3((nuclei.xmin+nuclei.xmax)/2,
						   (nuclei.ymin+nuclei.ymax)/2,
						   (nuclei.zmin+nuclei.zmax)/2);*/
		pm.setEulerAngles(4.26073, 1.18864, -3.54324);

		std::function<PointModel(vector<double>)> 
			unpackPoint = [=](vector<double> v) {
			PointModel m;
			m.ax      = v[0];
			m.rh      = v[1];
			
			m.metal.x = v[2];//4.165; //randomPointModel.metal.x;
			m.metal.y = v[3];//18.875; //randomPointModel.metal.y;
			m.metal.z = v[4];//17.180; //randomPointModel.metal.z;
						
			m.setEulerAngles(v[5],v[6],v[7]);
			
			return m;
		};

		std::function<vector<double>(PointModel)>
			packPoint = [](const PointModel& m){
			vector<double> vec;
			vec.push_back(m.ax);
			vec.push_back(m.rh);
			
			vec.push_back(m.metal.x);
			vec.push_back(m.metal.y);
			vec.push_back(m.metal.z);

			vec.push_back(m.angle_x);
			vec.push_back(m.angle_y);
			vec.push_back(m.angle_z);
			
			return vec;
		};

		std::function<PointModel(unsigned long) > randomPoint = [=](unsigned long seed) {
			PRNG prng(seed);
			RandomDist rand;

			PointModel m;
			m.ax = randomPointModel.ax; //(1-2*rand(prng));
			m.rh = randomPointModel.rh;//(1-2*rand(prng));

			m.metal.x = randomPointModel.metal.x;
			m.metal.y = randomPointModel.metal.y;
			m.metal.z = randomPointModel.metal.z;
			
			m.setEulerAngles(randomPointModel.angle_x,randomPointModel.angle_y,randomPointModel.angle_z);
			return m;
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
		gm.stddev = 0.2;

		std::function<GaussModel(vector<double>)>
			packGauss = [=](vector<double> v){
			GaussModel m;
			m.ax      = v[0];
			m.rh      = v[1];
			
			m.metal.x = v[2];//randomPointModel.metal.x;
			m.metal.y = v[3];//randomPointModel.metal.y;
			m.metal.z = v[4];//randomPointModel.metal.z;
												  
			m.setEulerAngles(v[5],v[6],v[7]);
												  
			m.stddev = v[8];
			return m;
		};

		std::function<vector<double>(const GaussModel&)>
			unpackGauss = [](const GaussModel& m){
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
