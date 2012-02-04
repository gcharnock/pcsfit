
/*
  Licenced under the GPL v2

  As an experiment, I decided to try to write this programe in a more
  C-ish style, even though its still C++ and a few C++ features are
  used where appropriate. I just tried to avoid using every feature in
  the C++ toolkit just because I could.
 */


#include <unistd.h>
#include <gsl/gsl_multimin.h>
#include <vector>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <sstream>
#include <cassert>
#include <boost/function.hpp>
#include <limits>
#include <fstream>

#include "fit.hpp"
#include "foreach.hpp"
#include "threads.hpp"
#include "model2.hpp"
#include "data.hpp"
#include "tests.hpp"

#define logMsg(x) flog << __FILE__ << "(" << __LINE__ << ")" << x << endl;

using namespace std;
using namespace boost;

bool main_on_iterate(const ErrorContext* context,unsigned long itN,gsl_multimin_fdfminimizer* min);


struct Options {
    Options ()
        : dont_rescale(false),
          dont_reg(false),
          integral_tol(1e-3),
          seed(0),
          scanparam(0),
          scanparam_min(0),
          scanparam_max(1),
          scanparam2(1),
          scanparam_min2(0),
          scanparam_max2(1),
          step_count(10),
          rel_error(1e-3) {
    }
    string params_file;
    string input_file;
    string method;

    string sketch3dFile;

    bool dont_rescale;
    bool dont_reg;
    double integral_tol;

    unsigned long seed;

    unsigned long scanparam;
    double scanparam_min;
    double scanparam_max;

    long freezeParam;

    long scanparam2;
    double scanparam_min2;
    double scanparam_max2;

    unsigned long step_count;

    double rel_error;
};

/********************************************************************************
 * Logging
 ********************************************************************************/

ofstream flog;


void run_scan(const Options& options,const ErrorContext* context,const ModelOptions* modelOptions,bool errorscan = false) {
    unsigned long size = context->model->size;
    double* params = (double*)alloca(size * sizeof(double));
    memcpy(params,context->params,size * sizeof(double));

    bool twoDScan = options.scanparam2 > -1;

    if(options.scanparam < 0 || options.scanparam > size) {
        cerr << "Invalid value for scanparam 1" << endl;
        return;
    }
    if(options.scanparam2 > long(size)) { //long(), shut up, gcc warrning
        cerr << "Invalid value for scanparam 2" << endl;
        return;
    }
    
    cout << "Paramiter set:" << endl;
    for(unsigned long i = 0; i < size; i++) { 
        cout << context->params[i];
        if(options.scanparam == i) {
            cout << "  - First scan dimension";
        }
        if(options.scanparam2 == long(i)) {
            cout << "  - Second scan dimension";
        }
        cout << endl;
    }
    cout << endl;

    for(unsigned long i = 0; i < options.step_count; i++) {
        double t1 = double(i)/options.step_count;
        double p1 = options.scanparam_min+(options.scanparam_max-options.scanparam_min)*t1;

        params[options.scanparam] = p1;

        unsigned long inner_count = twoDScan ? options.step_count : 1;
        for(unsigned long j = 0; j < inner_count; j++) {
            double p2 = 0; //Shut up the warning
            if(twoDScan) {
                double t2 = double(j)/options.step_count;
                p2 = options.scanparam_min2 + (options.scanparam_max2-options.scanparam_min2)*t2;
                params[options.scanparam2] = p2;
            }


            double value;
            double* gradient = (double*)alloca(size*sizeof(double));
            if(!errorscan) {
                context->model->modelf(Vec3d(10,50,0),params,&value,gradient,modelOptions);
            } else {
                eval_error(context,params,&value,gradient);
            }
            
            if(twoDScan) {
                cout << p1 << " " << p2 << " " << value;
            } else {
                cout << p1 << " " << value;
            }
            for(unsigned long i = 0; i < size; i++) {
                cout << " " << gradient[i];
            }
            cout << endl;
        }
    }
}

int main(int argc,char** argv) {
	cout.width(20);
	cout.fill(' ');
	cout.precision(17);
	cout << fixed << showpos;

	//Parse the command line and decide what to do
	using namespace boost::program_options;
	variables_map variablesMap;
    
    string command;
    Options options;
    options.seed = 0;

    bool rescale = true;
	try {
        positional_options_description posOpt;
        posOpt.add("command",-1);

		options_description optDesc("Options");

		optDesc.add_options()
            ("command",value<string>(&command),"fit|scan|errorscan|sketch3d|makedata|selftest")
            ("params-file,p",value<string>(&options.params_file),"")
			("input-file,i",value<string>(&options.input_file),"specify an input file")
			("help,h","print this help page and exit")
            ("method,m",value<string>(&options.method)->default_value("bfgs"),"")
            ("dont-reg","Don't regualize")
            ("freeze,f",value<long>(&options.freezeParam)->default_value(-1),"Don't vary a paramiter")
            ("dont-rescale","")
            ("sketch3d-file",value<string>(&options.sketch3dFile),"")
			("random-model-type",
			 "Specify the type of the model to use for generating random data")
			("seed,s",value<unsigned long>(&options.seed),"Specify a seed for the random number generator")
			("param-to-scan,-s",value<unsigned long>(&options.scanparam)->default_value(0),"")
			("param-to-scan2",value<long>(&options.scanparam2)->default_value(-1),"")
			("min",value<double>(&options.scanparam_min)->default_value(0),"")
			("max",value<double>(&options.scanparam_max)->default_value(1),"")
			("min2",value<double>(&options.scanparam_min2)->default_value(0),"")
			("max2",value<double>(&options.scanparam_max2)->default_value(1),"")
            ("steps",value<unsigned long>(&options.step_count)->default_value(20),"");
		try {
			store(command_line_parser(argc,argv).options(optDesc).positional(posOpt).run(),variablesMap);
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
        if(variablesMap.count("command")) {
            command = variablesMap["command"].as<string>();
        }
        if(!(command == "fit" || command == "errorscan" || command == "sketch3d" ||
             command == "scan" || command == "selftest" || command == "makedata")) {
			cout << "Usage:" << endl;
			cout << optDesc << endl;
			return -1;
        }


	} catch(std::exception& e) {
		//I'm not sure if this is possible
		cerr << "An error occured when trying to parse the command line" << endl;
	}
	notify(variablesMap);

    if(variablesMap.count("rescale") > 0) {
        rescale = false;
        options.dont_rescale = true;
    }

    if(command == "fit") {
        if(variablesMap.count("input-file") != 1 && 
           variablesMap.count("random-data") == 0) {
            cout << "Specify exactly one input file" << endl;
            return 0;
        }
    }


	//Open the log file and params.log file
	flog.open("log.log");
	if(!flog.is_open()) {
		cerr << "Could not open log.log for writing" << endl;
		return 1;
	}

	//Seed the random number generator
	PRNG prng;
	RandomDist dist;
    prng.seed(options.seed);
    logMsg("Random number generator seeded with " << options.seed);

    //Start the thread pool
    Multithreader<fdf_t> pool;


    //Set the model options, such as integral accuracy
    ModelOptions modelOptions;
    modelOptions.absError = 0;
    modelOptions.relError = 1e-2;


    //If we just want to run the self tests, do it and stop
    if(command == "selftest") {
        testMaths(prng);
        testModel(prng,&pool,&modelOptions);
        return 0;
    }
    //----------------------------------------------------------------------//
    // We're not self testing, so we need to decide of a model and get
    // paramiters for that model

	//Now, which model should we use?
	const Model* model;
    std::vector<double> params;

    double* params_start = NULL;
    double* params_opt   = NULL;


    if(command != "sketch3d") {
        int retVal =  parse_params_file(options.params_file,&model,&params);
        if(retVal != PARSE_SUCESS) {
            cerr << "Parse of param file failed: ";
            switch(retVal) {
            case UNKNOWN_MODEL:
                cerr << "unknown model" << endl;
                break;
            case NOT_ENOUGH_PARAMS:
                cerr << "not enough params" << endl;
                break;
            case PARAM_FILE_NOT_FOUND:
                cerr << "file " << options.params_file << " not found" << endl;
                break;
            default:
                assert(false);
            }
            return -1;
        }

    	//Set up buffers to store the paramiters
    	params_start         = (double*)alloca(model->size*sizeof(double));
    	params_opt           = (double*)alloca(model->size*sizeof(double));


    	for(unsigned long i = 0; i < model->size;i++) {
    		params_start[i] = params[i];
    	}
    }
    if(command == "scan") {
        ErrorContext context;
        context.params  = params_start;
        context.model   = model;
        context.pool    = &pool;
        run_scan(options,&context,false);
        return 0;
    }

    //----------------------------------------------------------------------//
    //We're not scanning as a paramiter varies, so we must be fitting,
    //so we need a dataset to do that.

	//Load the data
	Dataset dataset;

	//Okay, we don't want tests, we are doing a fitting. Where should
	//data dataset come from?

    if (command == "makedata") {
        random_data(prng,*model,params_start,20,&dataset,&modelOptions);
        for(unsigned long i = 0; i < dataset.nuclei.size(); i++) {
            cout << dataset.vals[i]       << " "
                 << dataset.nuclei[i].x() << " "
                 << dataset.nuclei[i].y() << " "
                 << dataset.nuclei[i].z() << endl;
        }
        return 0;
	} else {
        //A file
        logMsg("Loading data from file " << options.input_file);
		loadData(options.input_file,&dataset);
    }
    if (dataset.nuclei.size() == 0) {
        cout << "WARNING: dataset contains no spins" << endl;
    }


    //Do we want to produce a publishable quality visualisation?  Is
    //so, do it now and quit
    if(command == "sketch3d") {
        ofstream out(options.sketch3dFile.c_str());
        if(!out.is_open()) {
            cerr << "Could not open file " << options.sketch3dFile << endl;
            return 1;
        }
        toSketch3d(out,&dataset);
        return 0;
    }

    //Stuff our dataset, starting params and model in an ErrorContex
	ErrorContext context;
	context.dataset = &dataset;
	context.params  = params_start;
	context.model   = model;
	context.pool    = &pool;
    context.rescale = rescale;

    //Do we want to scan thoughZ
    if(command == "errorscan") {
        run_scan(options,&context,&modelOptions,true);
        return 0;
    }

	//Okay, we've got everything, do the fitting
	double errorFinal = 666; //Set these to easily recognisable uninitalised values

	do_fit_with_grad(&context,params_opt,&errorFinal,&main_on_iterate);

    for(unsigned long j = 0;j < model->size; j++) {
		cout << name_param(j) << " start =" << params_start[j]
			 << " final = " << params_start[j]*params_opt[j] << endl;
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
        double gradj = gsl_vector_get(g,j);
		norm += gradj*gradj;
        assert(isfinite(norm));
	}
    assert(isfinite(norm));

    cout << itN << "\t";
    for(unsigned long j = 0;j< x->size;j++) {cout << (gsl_vector_get(x,j)*context->params[j]) << "\t";}
	cout << "\tf(x) = " << fx << "\t|grad| = " << norm << endl;

	return itN < 6000 && norm > 1e-30;
}

