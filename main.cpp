#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>
#include <unistd.h>
#include <gsl/gsl_multimin.h>
#include <vector>


#include "model.hpp"
#include "data.hpp"


using namespace std;

/*********************************************************************************
 * Global State
 *********************************************************************************/

unsigned long length = 0;

//The current model
GaussModel gaussModel;

//Nuclear Positions
Nuclei nuclei;

//Experimental Vals
Vals expVals;

//Calculated values
Vals calcVals;

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



void* thread_main(void* threadIdPonter) {
    //Whoami? Am I the special root thread?
    long threadId = *(long*)threadIdPonter;


    barrierGuard(pthread_barrier_wait(&barr));
    while(!threadExit) {
        barrierGuard(pthread_barrier_wait(&barr));
        /* Master thread is setting the global state here */
        barrierGuard(pthread_barrier_wait(&barr));
        printf("thread %ld is working...\n",threadId);

        for(unsigned long i=0;i<length;i+=(numCPU)) {
			gaussModel.eval(nuclei[i].x,nuclei[i].y,nuclei[i].z,expVals[i]/100);
        }
    }
}

double minf(const gsl_vector * v, void *) {

    barrierGuard(pthread_barrier_wait(&barr));
    gaussModel.ax = gsl_vector_get(v, 0);
    gaussModel.rh = gsl_vector_get(v, 1);

    printf("%le, %le \n",gaussModel.ax,gaussModel.rh);
    /*Other threads do work here*/
    barrierGuard(pthread_barrier_wait(&barr));

    //Now reduce the results
    double total = 0;
	
    for(unsigned long i = 0;i<length;i++) {
        double diff = expVals[i] - calcVals[i];
        total += diff*diff;
    }
    printf("err = %le \n",total);
    printf("\n",total);

    return total;
}

int main() {
    //Load the data
	pairNucVals data = loadData("dataset_one.inp");
	
	nuclei = data.first;
	expVals = data.second;

	calcVals.resize(expVals.size());
	length = expVals.size();

	//Update the global state


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
