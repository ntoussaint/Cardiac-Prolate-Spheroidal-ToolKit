#include <vtkLandmarkWidget.h>

#include <vtkObjectFactory.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>

vtkStandardNewMacro(vtkLandmarkWidget);
vtkCxxRevisionMacro(vtkLandmarkWidget, "$Revision: 1315 $");


////----------------------------------------------------------------------------
//class SphereActorCallback : public vtkCommand
//{
//public:
//  static SphereActorCallback *New()
//  { return new SphereActorCallback; }

//  void Execute(vtkObject *caller,
//               unsigned long event,
//               void *vtkNotUsed(callData))
//  {
//    if (this->Actor != NULL)
//    {
//      vtkActor* caller = vtkActor::SafeDownCast (caller);
//      if (caller && (event == vtkCommand::ModifiedEvent))
//      {
//        this->Actor->SetVisibility(caller->GetVisibility());
//        this->Actor->SetMapper(caller->GetMapper());
//        this->Actor->Modified();
//      }
//    }
//  }

//  void SetActor (vtkActor* arg)
//  {
//    if (this->Actor == arg)
//      return;
//    if (this->Actor)
//      this->Actor->UnRegister (this);
//    this->Actor = arg;
//    if (this->Actor)
//      this->Actor->Register(this);
//  }

//  vtkActor* GetActor()
//  {
//    return this->Actor;
//  }

//protected:
//  SphereActorCallback()
//  {
//    this->Actor = NULL;
//  }
//  ~SphereActorCallback()
//  {
//    if (this->Actor)
//      this->Actor->Delete ();
//  }

//  vtkActor* Actor;

//};


vtkLandmarkWidget::vtkLandmarkWidget()
{
  this->Command = vtkLandmarkWidgetCommand::New();
  this->Command->SetLandmark (this);
  this->Value = 0.0;
}

vtkLandmarkWidget::~vtkLandmarkWidget()
{
  this->Command->Delete();
}

void vtkLandmarkWidgetCommand::Execute(vtkObject *   caller, 
				       unsigned long event, 
				       void *        callData)
{
  vtkSphereWidget* l = vtkSphereWidget::SafeDownCast(caller);
  
  if (event == vtkCommand::InteractionEvent)
  {
    if (this->Landmark)
    {
      this->Landmark->SetCenter (l->GetCenter());
      if (this->Landmark->GetInteractor())
	this->Landmark->GetInteractor()->Render();
    }
  }
  if ( (event == vtkCommand::EnableEvent) ||
       (event == vtkCommand::DisableEvent) )
  {
    if (this->Landmark)
      this->Landmark->SetEnabled (l->GetEnabled());
  }
}

void vtkLandmarkWidgetCommand::SetLandmark (vtkSphereWidget* l)
{
  this->Landmark = l;
}


void vtkLandmarkWidget::SetEnabled( int val)
{
  Superclass::SetEnabled( val);
//  vtkRenderWindowInteractor *i = this->Interactor;
//  if (!i || !i->GetRenderWindow())
//      return;
//  vtkRendererCollection* renderers = i->GetRenderWindow()->GetRenderers();
//  renderers->InitTraversal();
//  std::cout<<"starting..."<<std::endl;
//  while(vtkRenderer* r = renderers->GetNextItem())
//  {
//      SphereActorCallback* cbk = SphereActorCallback::New();
//      vtkActor* actor = vtkActor::New();
//      cbk->SetActor (actor);
//      actor->SetMapper (this->SphereActor->GetMapper());
//      actor->SetVisibility(this->SphereActor->GetVisibility());
//      actor->SetProperty(this->SphereActor->GetProperty());
//      r->AddActor(actor);
//      std::cout<<"adding actor "<<actor<<" to renderer : "<<r<<std::endl;
//      actor->Delete();
//  }
}
