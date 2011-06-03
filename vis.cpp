
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

using namespace std;

void* thread_main(void*) {
	long i = 0;
	while(true) {
		cout << "starting a long calculation" << endl;
		sleep(5);
		cout << "calculating..." << endl;
		sleep(5);
		cout << "Almost there..." << endl;
		sleep(5);
		i++;
		cout << "Done!" << endl;
	}
	return NULL;
}

class TimerCallback : public vtkCommand {
public:
	static TimerCallback* New() {
		return new TimerCallback;
	}
	virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId,void *vtkNotUsed(callData)) {
		if(vtkCommand::TimerEvent == eventId) {
			cout << "Timeout..." << endl;
		}
	}
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

	}

	~FittingWindow() {
	}
	void setNuclei(const Nuclei& nuclei) {
		mNuclei = nuclei;
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
		for(unsigned long i = 0; i<length;i++) {
			double c = (vals[i]-vals.min)/(vals.max-vals.min);

			vtkActor *sphereActor = vtkActor::New();
			sphereActor->SetMapper(mSphereMapper);
			sphereActor->SetScale(0.1*sqrt(abs(vals[i])));
			sphereActor->SetPosition(mNuclei[i].x,mNuclei[i].y,mNuclei[i].z);
			sphereActor->GetProperty()->SetColor(c,0,1-c);

			renderer->AddActor(sphereActor);
			
			sphereActor->Delete();
		}
	}

	//Threading stuff
    pthread_t thread;


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


GaussModel gaussModel;



int main () {

    //pair<Nuclei,Vals> pair_nv = fakeData(&gaussModel,100);
    pair<Nuclei,Vals> pair_nv = loadData("dataset_one.inp");
    Nuclei nuclei = pair_nv.first;
    Vals vals = pair_nv.second;
	Vals modelVals;
	modelVals.resize(vals.size());
	
	gaussModel.ax = 10000;
	gaussModel.metal.x = (nuclei.xmax + nuclei.xmin)/2;
	gaussModel.metal.y = (nuclei.ymax + nuclei.ymin)/2;
	gaussModel.metal.z = (nuclei.zmax + nuclei.zmin)/2;

	gaussModel.bulkEval(nuclei,modelVals);

	FittingWindow fittingWindow;
	fittingWindow.setNuclei(nuclei);
	fittingWindow.setCalcVals(modelVals);
	fittingWindow.setExpVals(vals);
	fittingWindow.setMetal(gaussModel.metal);
	fittingWindow.mainLoop();

	return 0;
}
