
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

	mAngle_x = 0;
	mAngle_y = 0;
	mAngle_z = 0;

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
    if(mDataChanged) {
        //There is new data avaiable. We should change the visualisation
        mDataChanged = false;
        update();
    }
}

void FittingWindow::setNuclei(const Nuclei& nuclei) {
    mNuclei = nuclei;
    updateNuclei(mCalcRenderer,mCalcSpheres);
    updateNuclei(mExpRenderer ,mExpSpheres);
}
void FittingWindow::setCalcVals(const Vals& calcVals,const Vector3& metal,
								double angle_x,double angle_y,double angle_z) {
    mCalcVals = calcVals;
    mMetal = metal;
    mDataChanged = true;

	mAngle_x = angle_x;
	mAngle_y = angle_y;
	mAngle_z = angle_z;

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
    updateVals(mCalcRenderer,mCalcVals,mCalcSpheres);
    updateVals(mExpRenderer, mExpVals ,mExpSpheres );

    //Move the arrow
    vtkActor* arrowActor = vtkActor::New();
    arrowActor->SetMapper(mArrowMapper);
    arrowActor->SetScale(4);
	arrowActor->RotateX(mAngle_x);
	arrowActor->RotateY(mAngle_y);
	arrowActor->RotateZ(mAngle_z);
    arrowActor->SetPosition(mMetal.x,mMetal.y,mMetal.z);
    arrowActor->GetProperty()->SetColor(1,1,0);
		
    mCalcRenderer->AddActor(arrowActor);

    arrowActor->Delete();

	mRenderWin->Render();
}

void FittingWindow::updateNuclei(vtkRenderer* renderer,
								 std::vector<vtkSmartPointer<vtkActor> >& actors) {
    cout << "FittingWindow::updateNuclei(" << renderer << "); mNuclei.size() = " << mNuclei.size() << endl;
    for(unsigned long i = 0; i<mNuclei.size();i++) {
        vtkSmartPointer<vtkActor> sphereActor = vtkActor::New();
        sphereActor->SetMapper(mSphereMapper);
        sphereActor->SetPosition(mNuclei[i].x,mNuclei[i].y,mNuclei[i].z);
        renderer->AddActor(sphereActor);
		actors.push_back(sphereActor);
    }
    renderer->ResetCameraClippingRange();
}

void FittingWindow::updateVals(vtkRenderer* renderer,
							   const Vals& vals,
							   std::vector<vtkSmartPointer<vtkActor> >& spheres) {
	for(unsigned long i = 0;i<vals.size();i++) {
		vtkSmartPointer<vtkActor> actor = spheres[i];
        double c = (vals[i]-mExpVals.min)/(mExpVals.max-mExpVals.min);
        actor->SetScale(0.1+0.1*sqrt(abs(vals[i])));
        actor->GetProperty()->SetColor(c,0,1-c);
	}
}

