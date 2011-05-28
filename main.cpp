#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>
#include <unistd.h>

#include "cuba.h"

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
long length = 0;
double* nuc_x = NULL;
double* nuc_y = NULL;
double* nuc_z = NULL;

double* val = NULL;

void loadData(const char* filename) {
	double x,y,z,v;

	FILE* fp = fopen(filename,"r");


	while(!feof(fp)) {
		if(!fscanf(fp,"%lg %lg %lg %lg \n",&x,&y,&z,&v)) {
            printf("Could not parse input\n");
            exit(1);
        }
		length++;
	}

	nuc_x = (double*)malloc(length * sizeof(double));
	nuc_y = (double*)malloc(length * sizeof(double));
	nuc_z = (double*)malloc(length * sizeof(double));

	val = (double*)malloc(length * sizeof(double));

	fseek(fp,0,SEEK_SET);

	long i = 0;
	while(!feof(fp)) {
		if(!fscanf(fp,"%lf %lg %lg %lg \n",&x,&y,&z,&v)) {
            printf("Could not parse input");
            exit(1);
        }
		nuc_x[i] = x;
		nuc_y[i] = y;
		nuc_z[i] = z;
		val[i] = v;
		i++;
	}
	assert(length == i);
	fclose(fp);
}

/*********************************************************************************
 * Global State
 *********************************************************************************/

//The dimensions of the cube we will be integrating over.
double cube_x_min, cube_x_max;
double cube_y_min, cube_y_max;
double cube_z_min, cube_z_max;

double xp;
double yp;
double zp;

double ax;
double rh;

double cutoff2;

double exponant;

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

/*********************************************************************************
 *  Integrand Function
 *********************************************************************************/


int Integrand(const int *ndim, const double xx[],
			  const int *ncomp, double ff[], void *userdata) {

	//Apply the transformation into the unit cube [0,1]^2
	double x = xx[0]*(cube_x_max - cube_x_min) + cube_x_min;
	double y = xx[1]*(cube_y_max - cube_y_min) + cube_y_min;
	double z = xx[2]*(cube_z_max - cube_z_min) + cube_z_min;

	double gx = xp - x;
	double gy = yp - y;
	double gz = zp - z;

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
	double f = (ax*(2*gz2 - gx2 - gy2) + rh*(3.0/2.0)*(gx2-gy2))/r5;

	double result = g * f;
	//Scale the result
	ff[0] = result * (cube_x_max - cube_x_min)*(cube_y_max - cube_y_min)*(cube_z_max - cube_z_min);

	return 0;
}

void* thread_main(void* threadIdPonter) {
    long threadId = *(long*)threadIdPonter;
    bool isRoot = threadId == 0;

    printf("Started thread %ld\n",threadId);

	timespec start;
	timespec end;

    barrierGuard(pthread_barrier_wait(&barr));
	if(isRoot) clock_gettime(CLOCK_MONOTONIC, &start);
    barrierGuard(pthread_barrier_wait(&barr));

	int comp, nregions, neval, fail;
	double integral[1], error[1], prob[1];

	for(long i=threadId; i< length;i+=(numCPU)) {
		xp = nuc_x[i];
		yp = nuc_y[i];
		zp = nuc_z[i];

		double epsAbs = val[i]/10;

        
		Cuhre(NDIM, NCOMP, Integrand, NULL,
			  EPSREL, epsAbs, VERBOSE | LAST,
			  MINEVAL, MAXEVAL, KEY,
			  &nregions, &neval, &fail, integral, error, prob);

		//printf("CUHRE RESULT:\tnregions  %d\tfail %d\n", nregions, fail);
        //printf("CUHRE RESULT:\t%e +- %e \t= %.3f\n",*integral, *error, *prob);

		//printf("nEval %d\n",neval);
	}
    barrierGuard(pthread_barrier_wait(&barr));
	if(isRoot) {
        clock_gettime(CLOCK_MONOTONIC, &end);
        // Calculate time it took
        double time_taken = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1e9;
        printf("Took %e\n",time_taken);
    }
	
	


}

int main() {

	loadData("dataset_one.inp");


	//Update the global state
	cube_x_min = -5; cube_x_max = 5;
	cube_y_min = -5; cube_y_max = 5;
	cube_z_min = -5; cube_z_max = 5;

	ax = 1;
	rh = 0; 

	cutoff2 = 0.05;

    // Barrier initialization

    numCPU = 1;
    //numCPU = sysconf(_SC_NPROCESSORS_ONLN);

    if(pthread_barrier_init(&barr, NULL, numCPU)) {
        printf("Could not create a barrier\n");
        return -1;
    }

    long* thread_args = (long*)malloc(numCPU * sizeof(long));
    pthread_t* threads = (pthread_t*)malloc(numCPU * sizeof(pthread_t));

    for(long i=0; i< numCPU; i++) {
        thread_args[i] = i;
        if(pthread_create(&threads[i], NULL, thread_main, (void *) &thread_args[i]) != 0) {
            printf("Could not create thread %ld\n", i);
            return -1;
        }
    }

    for(long i=0; i< numCPU; i++) {
        if(pthread_join(threads[i], NULL)) {
            printf("Could not join thread %ld\n", i);
            return -1;
        }
    }
    return 0;
}
