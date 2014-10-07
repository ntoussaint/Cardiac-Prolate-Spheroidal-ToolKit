#ifndef _vtk_LandmarkWidget_h_
#define _vtk_LandmarkWidget_h_

#include <vtkSphereWidget.h>
#include <vtkCommand.h>
#include <vtkSetGet.h>

//BTX
class VTK_EXPORT vtkLandmarkWidgetCommand : public vtkCommand
{    
 public:  
  static  vtkLandmarkWidgetCommand* New() { return new vtkLandmarkWidgetCommand; }
  void Execute(vtkObject *   caller, 
               unsigned long event, 
               void *        callData);
  void SetLandmark (vtkSphereWidget* l);

 protected:
  vtkLandmarkWidgetCommand()
  {
    this->Landmark = NULL;
  }
  ~vtkLandmarkWidgetCommand(){}  
  
 private:
  vtkSphereWidget* Landmark;
};

class VTK_EXPORT vtkLandmarkWidget : public vtkSphereWidget
{
 public:
  static vtkLandmarkWidget* New();
  vtkTypeRevisionMacro(vtkLandmarkWidget, vtkSphereWidget);
  vtkGetObjectMacro (Command, vtkLandmarkWidgetCommand);
  vtkGetMacro (Value, double);
  vtkSetMacro (Value, double);  

  vtkGetObjectMacro(SphereActor, vtkActor);
  vtkGetObjectMacro(HandleActor, vtkActor);

  virtual void SetEnabled(int);
 protected:
  vtkLandmarkWidget();
  ~vtkLandmarkWidget();
  
 private:
  vtkLandmarkWidgetCommand* Command;
  double Value;
  
};
//ETX

#endif //_vtk_LandmarkWidget_h_
