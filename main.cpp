#include <unistd.h>
#include <gsl/gsl_multimin.h>
#include <vector>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <sstream>
#include <cassert>
#include <boost/function.hpp>
#include <limits>

#include "fit.hpp"
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

bool main_on_iterate(const ErrorContext* context,unsigned long itN,gsl_multimin_fdfminimizer* min);

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

/********************************************************************************
 * Logging
 ********************************************************************************/

ofstream fout;
ofstream flog;


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
	unsigned long seed;
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
			("show-only","Show the data only, don't try and fit. Implies --gui")
			("seed",value<unsigned long>(&seed),"Specify a seed for the random number generator");

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

	//Open the log file and params.log file
	flog.open("log.log");
	if(!flog.is_open()) {
		cerr << "Could not open log.log for writing" << endl;
		return 1;
	}

	//Seed the random number generator
	PRNG prng;
	RandomDist dist;

	if(variablesMap.count("seed") > 1) {
		prng.seed(seed);
	}


    //Start the thread pool
    Multithreader<fdf_t> pool;

	//Load the data
	Dataset dataset;

	//If tests were requested, perform them and exit
	if(variablesMap.count("run-tests") > 0) {
		testModel(prng,&pool);
		return 0;
	}


	//Okay, we don't want tests, we are doing a fitting. Where should
	//data dataset come from?

	if(variablesMap.count("random-data") == 0) {
		logMsg("Opening the data");
		loadData(filename,&dataset);
	} else {
		//We want to generate it randomly. Which model should we use?
		//TODO
		logMsg("Generating a random point model");
		return -5;
    }

	//Now, which model should we use for fitting?
	const Model* model;

	if(modelType == "point") {
		model = &point_model;
	} else if(modelType == "gauss") {
		model = &gaussian_model;
	} else {
		cerr << "Unknown model type" << endl;
		return 1;
	}

	//What starting paramiters should we use?
	double* params_start       = (double*)alloca(model->size*sizeof(double));
	double* params_opt         = (double*)alloca(model->size*sizeof(double));

	//For now, let them be random
	for(unsigned long i = 0; i < model->size;i++) {
		params_start[i] = dist(prng);
	}
	params_start[PARAM_X] = 14.2815239342;
	params_start[PARAM_Y] = 58.2035463806;
	params_start[PARAM_Z] = 4.9245389075;

	params_start[PARAM_CHI1] = -10807.9638649971;
	params_start[PARAM_CHI2] = +118058.3614822003;
	params_start[PARAM_CHIXY] = +122257.0610884236;
	params_start[PARAM_CHIXZ] = +148620.2078244424;
	params_start[PARAM_CHIYZ] = +96579.0841891781;

	//But the stddev of the gaussian model should be small
	if(model == &gaussian_model) {params_start[PARAM_STDDEV] = 0.0000001;}

	//Open the params file TODO: What is the purpose of this file again?
	fout.open("params.log");
	if(!fout.is_open()) {
		cerr << "Could not open params.log for writing" << endl;
		return 1;
	}
	logMsg("Log file opened");

	//Okay, we've got everything, do the fitting

	ErrorContext context;
	context.dataset = &dataset;
	context.params  = params_start;
	context.model   = model;
	context.pool    = &pool;

	double errorFinal = 666; //Set these to easily recognisable uninitalised values

	do_fit_with_grad(&context,params_opt,&errorFinal,&main_on_iterate);
    
    for(unsigned long j = 0;j < model->size; j++) {
		cout << name_param(POINT_PARAM(j)) << " start =" << params_start[j]
			 << " final = " << params_opt[j] << endl;
	}
	cout.setf(ios::floatfield,ios::scientific);
	
	cout << "The final error was: " << errorFinal << endl;

	cout << "================================================================================" << endl;

    
    return 0;
}

//What to do on each itterations?

bool main_on_iterate(const ErrorContext* context,unsigned long itN,gsl_multimin_fdfminimizer* minimizer) {
	double fx = gsl_multimin_fdfminimizer_minimum(minimizer);
	gsl_vector* g = gsl_multimin_fdfminimizer_gradient(minimizer);
    gsl_vector* x = gsl_multimin_fdfminimizer_x(minimizer);

	double norm = 0;
	for(unsigned long j = 0; j < g->size; j++) {
		norm += gsl_vector_get(g,j)*gsl_vector_get(g,j);
	}

    cout << itN << "\t";
    for(unsigned long j = 0;j< x->size;j++) {cout << gsl_vector_get(x,j) << "\t";}
	cout << "\tf(x) = " << fx << "\t|grad| = " << norm << endl;

	return itN < 5000;
}

