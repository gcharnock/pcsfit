
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

class TimerCallback;

class FittingWindow {
    friend class TimerCallback;
public:
	FittingWindow();
    void onTimeout(unsigned long n);
	~FittingWindow();

    //These three functions set the data that the visualisation
    //visualised. As VTK is not thread safe the display doesn't update
    //imediately, instead these functions make a copy of the data
	void setNuclei(const Nuclei& nuclei);
	void setCalcVals(const Vals& calcVals,const Vector3& metal);
	void setExpVals(const Vals& expVals);

	void start();
	void mainLoop();
private:
    void update();
    void updateNuclei(vtkRenderer* renderer);
	void updateVals(vtkRenderer* renderer,const Vals& vals);
    void onTimeout();

	//Threading stuff
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

    bool    mDataChanged;
    Vals    mExpVals;
    Vals    mCalcVals;
	Vector3 mMetal;
	Nuclei  mNuclei;
	//Where to visualise the metal
};
