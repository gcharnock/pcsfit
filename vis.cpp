
#include <iostream>

#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

#include "data.hpp"
#include <math.h>

using namespace std;

int main () {
    //Data data = fakeData(100);
    pair<Nuclei,Vals> pair_nv = loadData("dataset_one.inp");
    Nuclei nuclei = pair_nv.first;
    Vals vals = pair_nv.second;

	// create sphere geometry
	vtkSphereSource *sphere = vtkSphereSource::New();
	sphere->SetRadius(1.0);
	sphere->SetCenter(0,0,0);
	sphere->SetThetaResolution(18);
	sphere->SetPhiResolution(18);

	// map to graphics library
	vtkPolyDataMapper *map = vtkPolyDataMapper::New();
	map->SetInput(sphere->GetOutput());

	vtkRenderer *renderer = vtkRenderer::New();

	// actor coordinates geometry, properties, transformation
    unsigned long length = nuclei.size();
    for(unsigned long i = 0; i<length;i++) {
        double c = (vals[i]-vals.min)/(vals.max-vals.min);
        cout << c << endl;

        vtkActor *sphereActor = vtkActor::New();
        sphereActor->SetMapper(map);
        sphereActor->SetScale(0.1*pow(abs(vals[i]),1.0/2));
        sphereActor->SetPosition(nuclei[i].x,nuclei[i].y,nuclei[i].z);
        sphereActor->GetProperty()->SetColor(c,0,1-c);

        renderer->AddActor(sphereActor);
    }

	// a render window

	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(renderer);

	// an interactor
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	// add the actor to the scene

	renderer->SetBackground(1,1,1); // Background color white

	// render an image (lights and cameras are created automatically)
	renWin->Render();

	// begin mouse interaction
	iren->Start();
  
	return 0;
}
