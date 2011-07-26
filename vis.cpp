
#include <iostream>


#include "unistd.h"

#include "data.hpp"
#include "model.hpp"
#include <math.h>

#include <gsl/gsl_multimin.h>

#include "vis.hpp"

#define PI 3.141592654

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
    mSphere->SetThetaResolution(5);
    mSphere->SetPhiResolution(5);

    mArrow = vtkArrowSource::New();

    //Creater the mappers
    mSphereMapper = vtkPolyDataMapper::New();
    mSphereMapper->SetInput(mSphere->GetOutput());

    mArrowMapper = vtkPolyDataMapper::New();
    mArrowMapper->SetInput(mArrow->GetOutput());

    mCalcRenderer = vtkRenderer::New();

    mRenderWin  = vtkRenderWindow::New();
    mRenderWin->AddRenderer(mCalcRenderer);

    mWindowInteractor = vtkRenderWindowInteractor::New();
    mRenderWin->SetInteractor(mWindowInteractor);
    mRenderWin->SetSize(800, 600);

    //We start by visualising iteraction 0
    visualisedIteraction = 0;

    mArrowActor = vtkActor::New();
    mArrowActor->SetMapper(mArrowMapper);
    mArrowActor->SetScale(4);
    mArrowActor->GetProperty()->SetColor(1,1,0);
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
    updateNuclei(mCalcSpheres,true);
    updateNuclei(mExpSpheres ,false);
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
    updateVals(mCalcRenderer,mExpVals ,mExpSpheres );

    //Move the arrow
    mArrowActor->SetOrientation(mAngle_x /2/PI * 360,
                                mAngle_y /2/PI * 360,
                                mAngle_z /2/PI * 360);
    mArrowActor->SetPosition(mMetal.x,mMetal.y,mMetal.z);
		
    mCalcRenderer->AddActor(mArrowActor);

	mRenderWin->Render();
}

void FittingWindow::updateNuclei(std::vector<vtkSmartPointer<vtkActor> >& actors,
                                 bool wireframe) {
    for(unsigned long i = 0; i<mNuclei.size();i++) {
        vtkSmartPointer<vtkActor> sphereActor = vtkActor::New();
        sphereActor->SetMapper(mSphereMapper);
        if(/*wireframe*/ i%2 == 0) {
            sphereActor->GetProperty()->SetRepresentationToWireframe();
        } else {
            sphereActor->GetProperty()->SetRepresentationToSurface();
        }


        sphereActor->SetPosition(mNuclei[i].x,mNuclei[i].y,mNuclei[i].z);
        mCalcRenderer->AddActor(sphereActor);
		actors.push_back(sphereActor);
    }
    mCalcRenderer->ResetCameraClippingRange();
}

void FittingWindow::updateVals(vtkRenderer* renderer,
							   const Vals& vals,
							   std::vector<vtkSmartPointer<vtkActor> >& spheres) {
    double maxAbs = abs(abs(mExpVals.max) > abs(mExpVals.min) ? mExpVals.max : mExpVals.min);
	for(unsigned long i = 0;i<vals.size();i++) {
		cout << vals[i] << endl;
		vtkSmartPointer<vtkActor> actor = spheres[i];

        double c = vals[i]/(maxAbs*2)+0.5;
        actor->SetScale(0.1+10*sqrt(abs(vals[i])/maxAbs));
        actor->GetProperty()->SetColor(0,c,1-c);
        if(/*wireframe*/ i%2 == 0) {
            actor->GetProperty()->SetRepresentationToWireframe();
        } else {
            actor->GetProperty()->SetRepresentationToSurface();
        }
	}
}

