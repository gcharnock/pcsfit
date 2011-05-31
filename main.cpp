#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>
#include <unistd.h>
#include <gsl/gsl_multimin.h>

#include <vector>

#include "cuba.h"

using namespace std;

#define NDIM 3
#define NCOMP 1
#define EPSREL 1
#define EPSABS 1e-8
#define VERBOSE 0
#define LAST 4
#define MINEVAL 0
#define MAXEVAL 50000

#define KEY 0

/********************************************************************************
 * Data
 ********************************************************************************/

struct DataPoint {
	double nuc_x;
	double nuc_y;
	double nuc_z;

	double val;
};

unsigned long length = 0;
vector<DataPoint> data;
vector<double> modelVals;

/*********************************************************************************
 * Global State
 *********************************************************************************/

struct Params {
	double xp;
	double yp;
	double zp;

	double ax;
	double rh;

	double exponant;
};

double cutoff2;

//The current model
Params params;

//The best model so far
Params bestParams;

//The dimensions of the cube we will be integrating over.
double cube_x_min, cube_x_max;
double cube_y_min, cube_y_max;
double cube_z_min, cube_z_max;


/********************************************************************************
 * Syncronisation
 ********************************************************************************/

// Barrier variable
pthread_barrier_t barr;

int numCPU;

void barrierGuard(int rc) {
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {
        printf("Could not wait on barrier\n");
        exit(-1);
    }
}

//When set to true, created threads should exit
bool threadExit = false;

/*********************************************************************************
 *  Load Data
 *********************************************************************************/

void loadData(const char* filename) {
	double x,y,z,v;

	FILE* fp = fopen(filename,"r");

    if(!fp) {
        printf("Could not open %s\n",filename);
        exit(1);
    }

	while(!feof(fp)) {
		if(!fscanf(fp,"%lf %lg %lg %lg \n",&x,&y,&z,&v)) {
            printf("Could not parse input");
            exit(1);
        }
		DataPoint p;
		p.nuc_x = x;
		p.nuc_y = y;
		p.nuc_z = z;
		p.val = v;
		data.push_back(p);
	}
	modelVals.resize(data.size());
	length = data.size();
	fclose(fp);
}

double rfloat() {
	return (double)rand()/(double)RAND_MAX;;
}

//Generates fake data by placing spins randomly in the [0,1]^3 cube
//and evaulating a random model. Good for testing how vunerable we are
//to local minima

void fakeData(unsigned long count) {
	srand(12345);

	for(unsigned long i=0;i<count;i++) {
		DataPoint p;
		p.nuc_x = rfloat();
		p.nuc_y = rfloat();
		p.nuc_z = rfloat();
		p.val = 0; /*TODO*/
		data.push_back(p);
	}
	length = data.size();
	modelVals.resize(data.size());
}


/*********************************************************************************
 *  Integrand Function
 *********************************************************************************/


int Integrand(const int *ndim, const double xx[],
			  const int *ncomp, double ff[], void *userdata) {

	double* nuclearLocation = (double*)userdata;

	//Apply the transformation into the unit cube [0,1]^2
	double x = xx[0]*(cube_x_max - cube_x_min) + cube_x_min;
	double y = xx[1]*(cube_y_max - cube_y_min) + cube_y_min;
	double z = xx[2]*(cube_z_max - cube_z_min) + cube_z_min;

	double gx = nuclearLocation[0] - x;
	double gy = nuclearLocation[1] - y;
	double gz = nuclearLocation[2] - z;

	double gx2 = gx*gx;
	double gy2 = gy*gy;
	double gz2 = gz*gz;

	double r2 = gx2+ gy2 + gz2;

	if(r2 < cutoff2) {
		ff[0] = 0;
		return 0;
	}

	double r = sqrt(r2);
	double r5 = r2*r2*r;

	double g = exp(-x*x-y*y-z*z);
	double f = (params.ax*(2*gz2 - gx2 - gy2) + params.rh*(3.0/2.0)*(gx2-gy2))/r5;

	double result = g * f;
	//Scale the result
	ff[0] = result * (cube_x_max - cube_x_min)*(cube_y_max - cube_y_min)*(cube_z_max - cube_z_min);

	return 0;
}

void* thread_main(void* threadIdPonter) {
    //Whoami? Am I the special root thread?
    long threadId = *(long*)threadIdPonter;

	int comp, nregions, neval, fail;
	double integral[1], error[1], prob[1];

    barrierGuard(pthread_barrier_wait(&barr));
    while(!threadExit) {
        barrierGuard(pthread_barrier_wait(&barr));
        /* Master thread is setting the global state here */
        barrierGuard(pthread_barrier_wait(&barr));
        printf("thread %ld is working...\n",threadId);

		unsigned long size = data.size();

        for(unsigned long i=0;i<size;i+=(numCPU)) {
			double nuclearLocation[3];
			nuclearLocation[0] = data[i].nuc_x;
            nuclearLocation[1] = data[i].nuc_y;
            nuclearLocation[2] = data[i].nuc_z;

            double epsAbs = data[i].val/1000;

        
            Cuhre(NDIM, NCOMP, Integrand, (void*)nuclearLocation,
                  EPSREL, epsAbs, VERBOSE | LAST,
                  MINEVAL, MAXEVAL, KEY,
                  &nregions, &neval, &fail, integral, error, prob);
            modelVals[i] = *integral;

            //printf("CUHRE RESULT:\tnregions  %d\tfail %d\n", nregions, fail);
            //printf("CUHRE RESULT:\t%e +- %e \t= %.3f\n",*integral, *error, *prob);

            //printf("nEval %d\n",neval);
        }
    }
}

double minf(const gsl_vector * v, void *) {

    barrierGuard(pthread_barrier_wait(&barr));
    params.ax = gsl_vector_get(v, 0);
    params.rh = gsl_vector_get(v, 1);

    printf("%le, %le \n",params.ax,params.rh);
    /*Other threads do work here*/
    barrierGuard(pthread_barrier_wait(&barr));

    //Now reduce the results
    double total = 0;
	
	unsigned long size = data.size();
	
    for(unsigned long i = 0;i<length;i++) {
        double diff = modelVals[i] - data[i].val;
        total += diff*diff;
    }
    printf("err = %le \n",total);
    printf("\n",total);

    return total;
}

int main() {
    //Load the data
	loadData("dataset_one.inp");


	//Update the global state
	cube_x_min = -5; cube_x_max = 5;
	cube_y_min = -5; cube_y_max = 5;
	cube_z_min = -5; cube_z_max = 5;

	params.ax = 1;
	params.rh = 0; 

	cutoff2 = 0.05;

    // pthreads initialization
    //numCPU = 1;
    numCPU = sysconf(_SC_NPROCESSORS_ONLN);
    if(pthread_barrier_init(&barr, NULL, numCPU + 1)) {
        printf("Could not create a barrier\n");
        return -1;
    }

    //Start the worker threads
    long* thread_args = (long*)malloc(numCPU * sizeof(long));
    pthread_t* threads = (pthread_t*)malloc(numCPU * sizeof(pthread_t));

    for(long i=0; i< numCPU; i++) {
        thread_args[i] = i;
        if(pthread_create(&threads[i], NULL, thread_main, (void *) &thread_args[i]) != 0) {
            printf("Could not create thread %ld\n", i);
            return -1;
        }
    }


    //Initalise the minimizer
    gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,2);
    gsl_multimin_function minfunc;

    minfunc.f = &minf;
    minfunc.n = 2;
    minfunc.params = NULL;

    gsl_vector *vec = gsl_vector_alloc (2);
    gsl_vector_set (vec, 0, 1.0);
    gsl_vector_set (vec, 1, 0.0);

    gsl_vector *step_size = gsl_vector_alloc (2);
    gsl_vector_set (step_size, 0, 0.1);
    gsl_vector_set (step_size, 1, 0.1);


    gsl_multimin_fminimizer_set (minimizer, &minfunc, vec, step_size);
     
    //When the treads are ready...

    barrierGuard(pthread_barrier_wait(&barr));
    for(long i = 0;i<20;i++) {
        //Root thread also handls timing
        timespec start;
        timespec end;
        clock_gettime(CLOCK_MONOTONIC, &start);


        gsl_multimin_fminimizer_iterate (minimizer);


        // Calculate time it took
        clock_gettime(CLOCK_MONOTONIC, &end);
        double time_taken = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1e9;
        //printf("Took %e\n",time_taken);
    }

    threadExit = true;

    return 0;
}
