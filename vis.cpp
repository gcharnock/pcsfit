
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"


int main () {
	// create sphere geometry
	vtkSphereSource *sphere = vtkSphereSource::New();
	sphere->SetRadius(1.0);
	sphere->SetCenter(0,0,0);
	sphere->SetThetaResolution(18);
	sphere->SetPhiResolution(18);

	// map to graphics library
	vtkPolyDataMapper *map = vtkPolyDataMapper::New();
	map->SetInput(sphere->GetOutput());
	//map->SetInput(sphere2->GetOutput());

	// actor coordinates geometry, properties, transformation
	vtkActor *aSphere = vtkActor::New();
	aSphere->SetMapper(map);
	aSphere->SetPosition(0,0,0);
	aSphere->GetProperty()->SetColor(1,0,0);

	vtkActor *aSphere2 = vtkActor::New();
	aSphere2->SetMapper(map);
	aSphere2->SetPosition(1,1,1);
	aSphere2->SetScale(0.3);
	aSphere2->GetProperty()->SetColor(0,0,1);


	// a renderer and render window
	vtkRenderer *ren1 = vtkRenderer::New();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(ren1);

	// an interactor
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	// add the actor to the scene
	ren1->AddActor(aSphere);
	ren1->AddActor(aSphere2);
	ren1->SetBackground(1,1,1); // Background color white

	// render an image (lights and cameras are created automatically)
	renWin->Render();

	// begin mouse interaction
	iren->Start();
  
	return 0;
}
