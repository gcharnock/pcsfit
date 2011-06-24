
#include <iostream>


#include "unistd.h"

#include "data.hpp"
#include "model.hpp"
#include <math.h>

#include <gsl/gsl_multimin.h>

#include "vis.hpp"

using namespace std;


/********************************************************************************
 * Visuals/main thread
 ********************************************************************************/

double iteraction_number = 0;

class FittingWindow;

class TimerCallback : public vtkCommand {
public:
	static TimerCallback* New() {
		return new TimerCallback;
	}
	virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId,void *vtkNotUsed(callData)) {
    if(vtkCommand::TimerEvent == eventId) {
        mParent->onTimeout(iteraction_number);
    }
}
    FittingWindow* mParent;
};


FittingWindow::FittingWindow() {
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
void FittingWindow::onTimeout(unsigned long n) {
    if(n != visualisedIteraction) {
        //There is new data avaiable. We should change the visualisation
        visualisedIteraction = n;
        mRenderWin->Render();
    }
}
FittingWindow::~FittingWindow() {
}
void FittingWindow::setNuclei(const Nuclei& nuclei) {
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
void FittingWindow::setCalcVals(const Vals& calcVals,const Vector3& metal) {
    mCalcVals = calcVals;
    setVals(mCalcRenderer, calcVals);

    mMetal = metal;
    vtkActor* arrowActor = vtkActor::New();
    arrowActor->SetMapper(mArrowMapper);
    arrowActor->SetScale(4);
    arrowActor->SetPosition(mMetal.x,mMetal.y,mMetal.z);
    arrowActor->GetProperty()->SetColor(1,1,0);
		
    mCalcRenderer->AddActor(arrowActor);

    arrowActor->Delete();

}
void FittingWindow::setExpVals(const Vals& expVals) {
    setVals(mExpRenderer, expVals);
}
void FittingWindow::start() {
    mWindowInteractor->Start();
}
void FittingWindow::mainLoop() {
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

void FittingWindow::setVals(vtkRenderer* renderer,const Vals& vals) {
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

