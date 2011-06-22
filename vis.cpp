
#include <iostream>

#include <vtkCommand.h>
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkArrowSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

#include "pthread.h"
#include "unistd.h"

#include "data.hpp"
#include "model.hpp"
#include <math.h>

#include <gsl/gsl_multimin.h>

using namespace std;


/*********************************************************************************
 * Global State
 *********************************************************************************/

//The current model
GaussModel gaussModel;

//Nuclear Positions
Nuclei nuclei;

//Experimental Vals
Vals expVals;

//Calculated values
volatile unsigned long iteraction_number = 0;
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

/********************************************************************************
 * Manager thread
 ********************************************************************************/

double minf(const gsl_vector * v, void *) {
    gaussModel.ax = gsl_vector_get(v, 0);
    gaussModel.rh = gsl_vector_get(v, 1);
    gaussModel.metal.x = gsl_vector_get(v, 2);
    gaussModel.metal.y = gsl_vector_get(v, 3);
    gaussModel.metal.z = gsl_vector_get(v, 4);
    gaussModel.setEulerAngles(gsl_vector_get(v, 5),gsl_vector_get(v, 6),gsl_vector_get(v, 7));
    gaussModel.exponant = gsl_vector_get(v, 8);

	gaussModel.bulkEval(nuclei,calcVals);

    cout << "x = " << gaussModel.metal.x;
    cout << " y = " << gaussModel.metal.y;
    cout << " z = " << gaussModel.metal.z;
    cout << " ax = " << gaussModel.ax;
    cout << " rh = " << gaussModel.rh;
    cout << " ax = " << gaussModel.angle_x;
    cout << " ay = " << gaussModel.angle_y;
    cout << " az = " << gaussModel.angle_z;
    cout << " exp = " << gaussModel.exponant;
    cout << endl;

    //Now reduce the results
    double total = 0;
	
    for(unsigned long i = 0;i<calcVals.size();i++) {
        double diff = expVals[i] - calcVals[i];
        total += diff*diff;
    }
    cout << "Error = " << total << endl;
    return total;
}


void* thread_main(void*) {
    numCPU = sysconf(_SC_NPROCESSORS_ONLN);

    static const long nParams = 9;

    //Initalise the minimizer
    gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,nParams);
    gsl_multimin_function minfunc;

    minfunc.f = &minf;
    minfunc.n = nParams;
    minfunc.params = NULL;

    gsl_vector *vec = gsl_vector_alloc (nParams);
    gsl_vector_set (vec, 0, gaussModel.ax); 
    gsl_vector_set (vec, 1, gaussModel.rh); 
    gsl_vector_set (vec, 2, gaussModel.metal.x);
    gsl_vector_set (vec, 3, gaussModel.metal.y);
    gsl_vector_set (vec, 4, gaussModel.metal.z);
    gsl_vector_set (vec, 5, gaussModel.angle_x);
    gsl_vector_set (vec, 6, gaussModel.angle_y);
    gsl_vector_set (vec, 7, gaussModel.angle_z);
    gsl_vector_set (vec, 8, gaussModel.exponant);


    gsl_vector *step_size = gsl_vector_alloc (nParams);
    gsl_vector_set (step_size, 0, 0.1);
    gsl_vector_set (step_size, 1, 0.1);
    gsl_vector_set (step_size, 2, 0.1);
    gsl_vector_set (step_size, 3, 0.1);
    gsl_vector_set (step_size, 4, 0.1);
    gsl_vector_set (step_size, 5, 0.1);
    gsl_vector_set (step_size, 6, 0.1);
    gsl_vector_set (step_size, 7, 0.1);
    gsl_vector_set (step_size, 8, 0.1);

    gsl_multimin_fminimizer_set (minimizer, &minfunc, vec, step_size);

    //When the treads are ready...

    for(long i = 0;true;i++) {
        //Root thread also handls timing
        timespec start;
        timespec end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        gsl_multimin_fminimizer_iterate (minimizer);
        iteraction_number = i;

        // Calculate time it took
        clock_gettime(CLOCK_MONOTONIC, &end);
        double time_taken = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1e9;
        printf("Took %e\n",time_taken);
    }
	return NULL;
}


/********************************************************************************
 * Visuals/main thread
 ********************************************************************************/

class FittingWindow;

class TimerCallback : public vtkCommand {
public:
	static TimerCallback* New() {
		return new TimerCallback;
	}
	virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId,void *vtkNotUsed(callData));
    FittingWindow* mParent;
};


class FittingWindow {
public:
	FittingWindow() {
		// create sphere geometry
		mSphere = vtkSphereSource::New();
		mSphere->SetRadius(1.0);
		mSphere->SetCenter(0,0,0);
		mSphere->SetThetaResolution(18);
		mSphere->SetPhiResolution(18);

		mArrow = vtkArrowSource::New();

		//Creater the mappers
		mSphereMapper = vtkPolyDataMapper::New();
		mSphereMapper->SetInput(mSphere->GetOutput());

		mArrowMapper = vtkPolyDataMapper::New();
		mArrowMapper->SetInput(mArrow->GetOutput());

		mCalcRenderer = vtkRenderer::New(); mCalcRenderer->SetViewport(0  ,0,0.5,1);
		mExpRenderer  = vtkRenderer::New(); mExpRenderer-> SetViewport(0.5,0,1  ,1);

		mRenderWin  = vtkRenderWindow::New();
		mRenderWin->AddRenderer(mCalcRenderer);
		mRenderWin->AddRenderer(mExpRenderer);

		mWindowInteractor = vtkRenderWindowInteractor::New();
		mRenderWin->SetInteractor(mWindowInteractor);
		mRenderWin->SetSize(800, 600);

        //We start by visualising iteraction 0
        visualisedIteraction = 0;


	}
    void onTimeout(unsigned long n) {
        if(n != visualisedIteraction) {
            //There is new data avaiable. We should change the visualisation
            visualisedIteraction = n;
            cout << "Need to update the visualisation" << endl;
            setCalcVals(calcVals);
            setMetal(gaussModel.metal);
            mRenderWin->Render();
        }
    }
	~FittingWindow() {
	}
	void setNuclei(const Nuclei& nuclei) {
		mNuclei = nuclei;
        unsigned long length = nuclei.size();
		for(unsigned long i = 0; i<length;i++) {
			vtkActor *sphereActor = vtkActor::New();
			sphereActor->SetMapper(mSphereMapper);
			sphereActor->SetPosition(mNuclei[i].x,mNuclei[i].y,mNuclei[i].z);
			mCalcRenderer->AddActor(sphereActor);
			
			sphereActor->Delete();
		}
		for(unsigned long i = 0; i<length;i++) {
			vtkActor *sphereActor = vtkActor::New();
			sphereActor->SetMapper(mSphereMapper);
			sphereActor->SetPosition(mNuclei[i].x,mNuclei[i].y,mNuclei[i].z);
			mExpRenderer->AddActor(sphereActor);
			
			sphereActor->Delete();
		}
	}
	void setCalcVals(const Vals& calcVals) {
		setVals(mCalcRenderer, calcVals);
	}
	void setExpVals(const Vals& expVals) {
		setVals(mExpRenderer, expVals);
	}
	void setMetal(const Vector3& metal) {
		mMetal = metal;
		vtkActor* arrowActor = vtkActor::New();
		arrowActor->SetMapper(mArrowMapper);
		arrowActor->SetScale(4);
		arrowActor->SetPosition(mMetal.x,mMetal.y,mMetal.z);
		arrowActor->GetProperty()->SetColor(1,1,0);
		
		mCalcRenderer->AddActor(arrowActor);

		arrowActor->Delete();
	}
	void start() {
		mWindowInteractor->Start();
	}
	void mainLoop() {
		mWindowInteractor->Initialize();

		vtkSmartPointer<TimerCallback> timerCallback = TimerCallback::New();
        timerCallback->mParent = this;
		mWindowInteractor->AddObserver(vtkCommand::TimerEvent,timerCallback);
		mWindowInteractor->CreateRepeatingTimer(500);

		pthread_t thread;

		if(pthread_create(&thread, NULL, thread_main, (void *)NULL) != 0) {
            printf("Could not create thread\n");
			return;
        }
		mWindowInteractor->Start();
	}
private:
	void setVals(vtkRenderer* renderer,const Vals& vals) {
		long length = vals.size();
        vtkActorCollection* actors = renderer->GetActors();
        vtkActor* actor = NULL;
        actors->InitTraversal();
        long i =0;
        while(true) {
            actor=actors->GetNextActor();
            if(!actor) {
                break;
            }
			double c = (vals[i]-vals.min)/(vals.max-vals.min);
			actor->SetScale(0.1*sqrt(abs(vals[i])));
			actor->GetProperty()->SetColor(c,0,1-c);
            i++;
        }

	}

	//Threading stuff
    pthread_t thread;
    unsigned long visualisedIteraction;


	//VTK Stuff
	vtkSmartPointer<vtkSphereSource> mSphere;
	vtkSmartPointer<vtkArrowSource> mArrow;

	vtkSmartPointer<vtkPolyDataMapper> mSphereMapper;
	vtkSmartPointer<vtkPolyDataMapper> mArrowMapper;

	vtkSmartPointer<vtkRenderer> mCalcRenderer;
	vtkSmartPointer<vtkRenderer> mExpRenderer;
	vtkSmartPointer<vtkRenderWindow> mRenderWin;
	vtkSmartPointer<vtkRenderWindowInteractor> mWindowInteractor;

	Nuclei mNuclei;
	//Where to visualise the metal
	Vector3 mMetal;

};
void TimerCallback::Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId,void *vtkNotUsed(callData)) {
    if(vtkCommand::TimerEvent == eventId) {
        mParent->onTimeout(iteraction_number);
    }
}



int main () {

    //pair<Nuclei,Vals> pair_nv = fakeData(&gaussModel,100);
    pair<Nuclei,Vals> pair_nv = loadData("dataset_one.inp");
    nuclei = pair_nv.first;
    expVals = pair_nv.second;
	calcVals.resize(expVals.size());
	
	gaussModel.ax = 0;
	gaussModel.metal.x = (nuclei.xmax + nuclei.xmin)/2;
	gaussModel.metal.y = (nuclei.ymax + nuclei.ymin)/2;
	gaussModel.metal.z = (nuclei.zmax + nuclei.zmin)/2;

    double sizex = nuclei.xmax - nuclei.xmin;
    double sizey = nuclei.ymax - nuclei.ymin;
    double sizez = nuclei.zmax - nuclei.zmin;

	gaussModel.cube_x_min = -sizex*2; gaussModel.cube_x_max = sizex*2;
	gaussModel.cube_y_min = -sizey*2; gaussModel.cube_y_max = sizey*2;
	gaussModel.cube_z_min = -sizez*2; gaussModel.cube_z_max = sizez*2;


    gaussModel.setEulerAngles(0,0,0);

	gaussModel.bulkEval(nuclei,calcVals);

	FittingWindow fittingWindow;
	fittingWindow.setNuclei(nuclei);
	fittingWindow.setCalcVals(calcVals);
	fittingWindow.setExpVals(expVals);
	fittingWindow.setMetal(gaussModel.metal);
	fittingWindow.mainLoop();

	return 0;
}
