
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


class TimerCallback : public vtkCommand {
public:
	static TimerCallback* New() {
		return new TimerCallback;
	}
	virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId,void *vtkNotUsed(callData)) {
        if(vtkCommand::TimerEvent == eventId) {
            mParent->onTimeout();
        }
    }
    FittingWindow* mParent;
};


FittingWindow::FittingWindow()
    : mDataChanged(true) {
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

FittingWindow::~FittingWindow() {
}

void FittingWindow::onTimeout() {
    cout << "Timeout" << endl;
    cout << "mDataChanged=" << mDataChanged << endl;
    if(mDataChanged) {
        cout << "we need to update" << endl;
        //There is new data avaiable. We should change the visualisation
        mDataChanged = false;
        update();
    }
}

void FittingWindow::setNuclei(const Nuclei& nuclei) {
    mNuclei = nuclei;
    updateNuclei(mCalcRenderer);
    updateNuclei(mExpRenderer);
}
void FittingWindow::setCalcVals(const Vals& calcVals,const Vector3& metal) {
    cout << "calc vals set, calcVals.size() = " << calcVals.size() << endl;
    mCalcVals = calcVals;
    mMetal = metal;
    mDataChanged = true;
    cout << "mDataChanged=" << mDataChanged << endl;    
}
void FittingWindow::setExpVals(const Vals& expVals) {
    mExpVals = expVals;
    mDataChanged = true;
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

    mWindowInteractor->Start();
}

void FittingWindow::update() {
    cout << "FittingWindow::update()" << endl;
    updateVals(mCalcRenderer,mCalcVals);
    updateVals(mExpRenderer, mExpVals);

    //Move the arrow
    vtkActor* arrowActor = vtkActor::New();
    arrowActor->SetMapper(mArrowMapper);
    arrowActor->SetScale(4);
    arrowActor->SetPosition(mMetal.x,mMetal.y,mMetal.z);
    arrowActor->GetProperty()->SetColor(1,1,0);
		
    mCalcRenderer->AddActor(arrowActor);

    arrowActor->Delete();
}

void FittingWindow::updateNuclei(vtkRenderer* renderer) {
    cout << "FittingWindow::updateNuclei(" << renderer << "); mNuclei.size() = " << mNuclei.size() << endl;
    for(unsigned long i = 0; i<mNuclei.size();i++) {
        vtkSmartPointer<vtkActor> sphereActor = vtkActor::New();
        sphereActor->SetMapper(mSphereMapper);
        sphereActor->SetPosition(mNuclei[i].x,mNuclei[i].y,mNuclei[i].z);
        renderer->AddActor(sphereActor);
    }
    renderer->ResetCameraClippingRange();
}

void FittingWindow::updateVals(vtkRenderer* renderer,const Vals& vals) {
    cout << "FittingWindow::updateVals(" << renderer << ",vals.size = " << vals.size() << ")" << endl;
    vtkActorCollection* actors = renderer->GetActors();
    vtkActor* actor = NULL;
    actors->InitTraversal();
    unsigned long i =0;
    while(true) {
        actor=actors->GetNextActor();
        if(!actor || vals.size() >= i) {
            break;
        }
        double c = (vals[i]-vals.min)/(vals.max-vals.min);
        actor->SetScale(0.1*sqrt(abs(vals[i])));
        cout << "Setting scale to " << 0.1*sqrt(abs(vals[i])) << endl;
        actor->GetProperty()->SetColor(c,0,1-c);
        i++;
    }
}

