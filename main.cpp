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
#include <sstream>

#include "foreach.hpp"
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
 * Logging
 ********************************************************************************/

ofstream fout;

/********************************************************************************
 * minimiser class
 ********************************************************************************/

template<typename T>
class Minimiser {
public:
	typedef std::vector<double> PList;
	Minimiser(boost::function<T(const PList&)> unpack,
			  boost::function<PList(const T&)> pack,  
			  boost::function<double(T)> min_funcion,		  
			  boost::function<void(const T&)> logger)		  
		: mUnpack(unpack),mPack(pack),mFMin(min_funcion),mLogger(logger) {
	}

	void minimise(T startingModel) {
		//Unpack the starting model
		PList plist = mPack(startingModel);
		unsigned long nParams = plist.size();

		//Initalise the minimiser
		gsl_multimin_fminimizer* gslmin =
			gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,nParams);

		gsl_multimin_function minfunc;
		minfunc.f = &Minimiser::minf_static;
		minfunc.n = nParams;
		minfunc.params = (void*)this;

		gsl_vector* vec = pList2GSLVec(plist);
		gsl_vector* step_size = gsl_vector_alloc(plist.size());
		
		for(unsigned long i = 0;i<plist.size();i++) {
			gsl_vector_set(step_size,i,plist[i] == 0 ? 1 : plist[i]*0.5);
		}

		gsl_multimin_fminimizer_set (gslmin, &minfunc, vec, step_size);
		for(unsigned long i = 0; i<200;i++) {
			gsl_multimin_fminimizer_iterate (gslmin);
		}

		gsl_multimin_fminimizer_free(gslmin);
		gsl_vector_free(step_size);
	}
private:
	double _minf(const gsl_vector * v) {
		T params = mUnpack(gSLVec2pList(v));
		mLogger(params);
		return mFMin(params);
	}
	static double minf_static(const gsl_vector * v, void* _this) {
		return ((Minimiser*)_this)->_minf(v);
	}

	static gsl_vector* pList2GSLVec(PList pList) {
		gsl_vector* vec = gsl_vector_alloc (pList.size());
		for(unsigned long i = 0; i<pList.size();i++) {
			gsl_vector_set(vec,i,pList[i]);
		}
		return vec;
	}
	static PList gSLVec2pList(const gsl_vector* vec) {
		PList retVal;
		for(size_t i = 0;i<vec->size;i++) {
			retVal.push_back(gsl_vector_get(vec,i));
		}
		return retVal;
	}

	boost::function<T(const PList&)> mUnpack;
	boost::function<PList(const T&)> mPack;
	boost::function<double(T)> mFMin;
	boost::function<void(T)> mLogger;	  
};

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
		 << thisModel.exponant << endl;

	stringstream fnameStream;
	fnameStream << "results/points" << models.size() << ".log";
	ofstream fpoints(fnameStream.str().c_str());
	for(unsigned long i = 0;i<results.size();i++) {
		fpoints << expVals[i] << " " << results[i] << " " << (expVals[i]-results[i]) << endl;
	}

    return total;
}

void onIterate(const GaussModel& m) {
	cout << "logging functions go here..." << endl;
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
	Minimiser<GaussModel> minimiser(&unpackGaussModel,
									&packGaussModel,
									&minf,
									&onIterate);
	
	//Set up the model
	GaussModel gm;
	gm.ax = 100.0;
	gm.rh = 0.0;
	gm.metal = Vector3((nuclei.xmin+nuclei.xmax)/2,
					   (nuclei.ymin+nuclei.ymax)/2,
					   (nuclei.zmin+nuclei.zmax)/2);
	gm.setEulerAngles(0,0,0);
	gm.exponant = 1;

    //Start the visualisation thread
    VisualThread visualThread;
    visualThread.fw.setNuclei(nuclei);
    visualThread.fw.setExpVals(expVals);
    boost::thread boostVisualThread(boost::ref(visualThread));

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
