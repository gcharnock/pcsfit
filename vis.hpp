
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
	void setCalcVals(const Vals& calcVals,const Vector3& metal,
					 double angle_x,double angle_y,double angle_z);
	void setExpVals(const Vals& expVals);

	void start();
	void mainLoop();
private:
    void update();
    void updateNuclei(std::vector<vtkSmartPointer<vtkActor> >& spheres,
                      bool wireframe);
	void updateVals(vtkRenderer* renderer,
					const Vals& vals,
					std::vector<vtkSmartPointer<vtkActor> >& spheres);
    void onTimeout();

	//Threading stuff
    unsigned long visualisedIteraction;


	//VTK Stuff
	vtkSmartPointer<vtkSphereSource> mSphere;
	vtkSmartPointer<vtkArrowSource> mArrow;

	vtkSmartPointer<vtkPolyDataMapper> mSphereMapper;
	vtkSmartPointer<vtkPolyDataMapper> mArrowMapper;

	vtkSmartPointer<vtkRenderer> mCalcRenderer;
	vtkSmartPointer<vtkRenderWindow> mRenderWin;
	vtkSmartPointer<vtkRenderWindowInteractor> mWindowInteractor;

	//VTK iteration over collections is pretty awful. We'll store
	//references to the spheres we create in std::vectors
	std::vector<vtkSmartPointer<vtkActor> > mCalcSpheres;
	std::vector<vtkSmartPointer<vtkActor> > mExpSpheres ;

    vtkSmartPointer<vtkActor> mArrowActor;

    bool    mDataChanged;
    Vals    mExpVals;
    Vals    mCalcVals;
	double  mAngle_x,mAngle_y,mAngle_z;
	Vector3 mMetal;
	Nuclei  mNuclei;
	//Where to visualise the metal
};
