#include <cmath>

#include "vtkImageView.h"
#include "vtkObjectFactory.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkMatrixToLinearTransform.h"

#include "vtkImageData.h"
#include "vtkPointSet.h"
#include <vtkInformation.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyle.h"
#include "vtkImageViewCornerAnnotation.h"
#include "vtkTextProperty.h"
#include "vtkCamera.h"
#include "vtkDataSetCollection.h"
#include "vtkProp3DCollection.h"

#include "vtkLookupTable.h"
#include "vtkColorTransferFunction.h"
#include "vtkPiecewiseFunction.h"
#include "vtkImageMapToColors.h"
#include "vtkScalarBarActor.h"
#include <vtkImageReslice.h>
#include <vtkSmartPointer.h>

#include "vtkCommand.h"
#include "vtkImageView2DCommand.h"

vtkCxxRevisionMacro(vtkImageView, "$Revision: 1 $");
// vtkStandardNewMacro(vtkImageView); // pure virtual class

vtkImageView::vtkImageView()
{
  this->OrientationMatrix       = vtkMatrix4x4::New();
  this->InvertOrientationMatrix = vtkMatrix4x4::New();
  this->OrientationTransform    = vtkMatrixToLinearTransform::New();
  this->CornerAnnotation        = vtkImageViewCornerAnnotation::New();
  this->TextProperty            = vtkTextProperty::New();
  this->DataSetCollection       = vtkDataSetCollection::New();
  this->DataSetActorCollection  = vtkProp3DCollection::New();
  this->LookupTable             = vtkLookupTable::New();
  this->WindowLevel             = vtkImageMapToColors::New();
  this->ScalarBar               = vtkScalarBarActor::New();
  
  this->Renderer               = 0;
  this->RenderWindow           = 0;
  this->Interactor             = 0;
  this->Input                  = 0;
  this->InternalMTime          = 0;
  this->InteractorStyle        = 0;
  this->IsInteractorInstalled  = 0;
  
  this->CornerAnnotation->SetNonlinearFontScaleFactor (0.3);
  this->CornerAnnotation->SetTextProperty ( this->TextProperty );
  this->CornerAnnotation->SetMaximumFontSize (46);
  this->CornerAnnotation->SetImageView (this);
  this->CornerAnnotation->PickableOff();
  
  this->ScalarBar->SetLabelFormat ("%.0f");
  this->ScalarBar->SetNumberOfLabels (3);
  this->ScalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  this->ScalarBar->SetLabelTextProperty (this->TextProperty);
  this->ScalarBar->GetLabelTextProperty()->SetFontSize (1);
  this->ScalarBar->GetLabelTextProperty()->BoldOff();
  this->ScalarBar->GetLabelTextProperty()->ShadowOff();
  this->ScalarBar->GetLabelTextProperty()->ItalicOff();
  this->ScalarBar->SetWidth (0.1);
  this->ScalarBar->SetHeight (0.5);
  this->ScalarBar->SetPosition (0.9,0.3);
  this->ScalarBar->PickableOff();
  this->ScalarBar->VisibilityOn();
  
  for(int i=0; i<3; i++)
    this->CurrentPoint[i] = 0.0; //VTK_DOUBLE_MIN;
  
  this->ShowAnnotations = true;
  this->ShowScalarBar = true;
  
  this->OrientationTransform->SetInput (this->OrientationMatrix);
  
  this->ColorWindow = 1;
  this->ColorLevel  = 0.5;

  this->WindowLevel->SetLookupTable( this->LookupTable );  
  this->WindowLevel->SetOutputFormatToRGB();
  this->ScalarBar->SetLookupTable( this->WindowLevel->GetLookupTable() );

  this->LookupTable->SetTableRange (0, 1);
  this->LookupTable->SetSaturationRange (0, 0);
  this->LookupTable->SetHueRange (0, 0);
  this->LookupTable->SetValueRange (0, 1);
  this->LookupTable->SetAlphaRange (0, 1);
  this->LookupTable->Build();
}

vtkImageView::~vtkImageView()
{
  this->OrientationTransform->SetInput ( NULL );
  
  this->OrientationMatrix->Delete();  
  this->InvertOrientationMatrix->Delete();  
  this->OrientationTransform->Delete();  
  
  this->CornerAnnotation->Delete();
  this->TextProperty->Delete();

  this->DataSetCollection->Delete();
  this->DataSetActorCollection->Delete();
  
  this->LookupTable->Delete();
  
  this->ScalarBar->Delete();
  this->WindowLevel->Delete();
  
  if( this->RenderWindow )
  {
    this->RenderWindow->Delete();
    this->RenderWindow = 0;
  }
  if( this->Renderer )
  {
      this->Renderer->Delete();
      this->Renderer = 0;
  }
  if( this->Interactor )
  {
    this->Interactor->Delete();
    this->Interactor = 0;
  }
  if (this->InteractorStyle)
  {
    this->InteractorStyle->Delete();
    this->InteractorStyle = 0;
  }
  if (this->Input)
  {
    this->Input->Delete();
    this->Input = 0;
  }
  
  std::cout<<"deleting a view. done"<<std::endl;
  
}

//----------------------------------------------------------------------------
unsigned long vtkImageView::GetMTime()
{
    typedef unsigned long MTimeType;

    MTimeType mTime = Superclass::GetMTime();

    vtkObject * objectsToInclude[] = {
        Renderer, RenderWindow, Interactor,
        InteractorStyle, WindowLevel, OrientationTransform, ScalarBar, OrientationMatrix,
        InvertOrientationMatrix, CornerAnnotation, TextProperty, LookupTable,
        ScalarBar, Input };

        const int numObjects = sizeof(objectsToInclude) / sizeof(vtkObject *);

        for ( int i(0); i<numObjects; ++i ) {
            vtkObject * object = objectsToInclude[i];
            if (object) {
                const MTimeType testMtime = object->GetMTime();
                if ( testMtime > mTime )
                    mTime = testMtime;
            }
        }
        return mTime;
}

//----------------------------------------------------------------------------
void vtkImageView::SetupInteractor(vtkRenderWindowInteractor *arg)
{
  this->UnInstallPipeline();
  
  vtkSetObjectBodyMacro (Interactor, vtkRenderWindowInteractor, arg);
  
  this->InstallPipeline();
}

//----------------------------------------------------------------------------
void vtkImageView::SetRenderWindow(vtkRenderWindow *arg)
{
  this->UnInstallPipeline();
  
  vtkSetObjectBodyMacro (RenderWindow, vtkRenderWindow, arg);
  
  if (this->RenderWindow && this->RenderWindow->GetInteractor())
  {
    this->SetupInteractor (this->RenderWindow->GetInteractor());
  }  
  this->InstallPipeline();
}

//----------------------------------------------------------------------------
void vtkImageView::SetRenderer(vtkRenderer *arg)
{
  this->UnInstallPipeline();
  
  vtkSetObjectBodyMacro (Renderer, vtkRenderer, arg);
  
  this->InstallPipeline();
}

//----------------------------------------------------------------------------
void vtkImageView::Render()
{
  if (this->RenderWindow)
  {
    if (!this->RenderWindow->GetNeverRendered())
    {
      if( this->GetMTime()>this->InternalMTime )
      {
	this->RenderWindow->Render();
        this->InternalMTime = this->GetMTime();
      }
    }
    else
    {
      this->RenderWindow->Render();
    }
  }
}

//----------------------------------------------------------------------------
void vtkImageView::SetInput(vtkImageData *arg, vtkMatrix4x4 *matrix) 
{  
  vtkSetObjectBodyMacro (Input, vtkImageData, arg);
  this->WindowLevel->SetInputData(arg);
  if (matrix) this->SetOrientationMatrix(matrix);
}

//----------------------------------------------------------------------------
void vtkImageView::SetInputConnection(vtkAlgorithmOutput* arg, vtkMatrix4x4 *matrix) 
{
  this->WindowLevel->SetInputConnection(arg);
  if (matrix) this->SetOrientationMatrix(matrix);
}

//----------------------------------------------------------------------------
void vtkImageView::InstallPipeline()
{
  if (this->RenderWindow && this->Renderer && !this->RenderWindow->HasRenderer (this->Renderer))
  {
    this->RenderWindow->AddRenderer(this->Renderer);
  }
  
  if (this->Interactor)
  {
    this->Interactor->SetRenderWindow(this->RenderWindow);
  }
  
  if (this->Renderer)
  {
    this->Renderer->AddViewProp ( this->CornerAnnotation );
    this->Renderer->AddViewProp ( this->ScalarBar );
    
    this->Renderer->GetActiveCamera()->ParallelProjectionOn();
  }
  
  this->InstallInteractor();
}

//----------------------------------------------------------------------------
void vtkImageView::UnInstallPipeline()
{
  this->UnInstallInteractor();
  
  if (this->Renderer)
  {
    this->Renderer->RemoveViewProp(this->CornerAnnotation);
    this->Renderer->RemoveViewProp(this->ScalarBar);
  }
  
  if (this->RenderWindow && this->Renderer)
  {
    this->RenderWindow->RemoveRenderer(this->Renderer);
  }
}

//----------------------------------------------------------------------------
void vtkImageView::GetWithinBoundsPosition (double* pos1, double* pos2)
{
  for (unsigned int i=0; i<3; i++) pos2[i] = pos1[i];
  
  if (!this->GetInput())
    return;
  
  int indices[3];
  this->GetImageCoordinatesFromWorldCoordinates (pos1, indices);
  
  this->GetInputAlgorithm()->UpdateInformation();
  int* w_extent = this->GetInputAlgorithm()->GetOutputInformation(0)->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
  bool out_of_bounds = false;
  
  for (unsigned int i=0; i<3; i++)
  {
    if (indices[i] > w_extent[2 * i + 1])
    {
      indices[i] = w_extent[2 * i + 1];
      out_of_bounds=true;
    }
    if (indices[i] < w_extent[2 * i])
    {
      indices[i] = w_extent[2 * i];
      out_of_bounds=true;
    }
  }
  
  if (out_of_bounds)
    this->GetWorldCoordinatesFromImageCoordinates (indices, pos2);
  else
    pos2 = pos1;
}


//----------------------------------------------------------------------------
void vtkImageView::SetCurrentPoint (double pos[3])
{
  double inside_pos[3];
  this->GetWithinBoundsPosition (pos, inside_pos);
  
  this->CurrentPoint[0] = inside_pos[0];
  this->CurrentPoint[1] = inside_pos[1];
  this->CurrentPoint[2] = inside_pos[2];
  this->InvokeEvent (vtkImageView::CurrentPointChangedEvent, NULL);
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkImageView::ResetCurrentPoint (void)
{
  if (!this->GetInput())
    return;
  
  int *wholeExtent = this->GetInput()->GetExtent();
  
  int center[3];
  for (unsigned int i=0; i<3; i++)
    center[i] = vtkMath::Round((double)(wholeExtent [2*i+1] - wholeExtent[2*i])/2.0);
  double position[3];
  this->GetWorldCoordinatesFromImageCoordinates (center, position);
  this->SetCurrentPoint (position);
}

//----------------------------------------------------------------------------
void vtkImageView::SetOrientationMatrix (vtkMatrix4x4* matrix)
{
  if (!matrix)
    return;
  
  vtkMatrix4x4* matrixcopy = vtkMatrix4x4::New();
  matrixcopy->DeepCopy (matrix);
  vtkSetObjectBodyMacro (OrientationMatrix, vtkMatrix4x4, matrixcopy);
  this->OrientationTransform->SetInput (this->OrientationMatrix);
  vtkMatrix4x4::Invert (this->OrientationMatrix, this->InvertOrientationMatrix);
  
  matrixcopy->Delete();
}

//----------------------------------------------------------------------------
void vtkImageView::SetLookupTable (vtkLookupTable* lookuptable)
{
  vtkSetObjectBodyMacro (LookupTable, vtkLookupTable, lookuptable);
  
  this->WindowLevel->SetLookupTable( this->LookupTable );
  this->ScalarBar->SetLookupTable( this->LookupTable );
}
//----------------------------------------------------------------------------
void vtkImageView::SetColorWindow(double s) 
{
  if (s < 0)
    s = 0;
  
  if (s == this->ColorWindow)
    return;
  
  this->ColorWindow = s;

  double range[2];
  range[0] = this->ColorLevel - this->ColorWindow / 2.0;
  range[1] = this->ColorLevel + this->ColorWindow / 2.0;  
  this->LookupTable->SetRange (range);
  this->LookupTable->SetRange (range);
  
  this->InvokeEvent (vtkImageView::WindowLevelChangedEvent, NULL);
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkImageView::SetColorLevel(double s) 
{
  if (s == this->ColorLevel)
    return;
  
  this->ColorLevel = s;

  double range[2];
  range[0] = this->ColorLevel - this->ColorWindow / 2.0;
  range[1] = this->ColorLevel + this->ColorWindow / 2.0;  
  this->LookupTable->SetRange (range);

  this->InvokeEvent (vtkImageView::WindowLevelChangedEvent, NULL);
  this->Modified();
}

void vtkImageView::SetColorRange( double r[2] )
{
  double level  = 0.5 * ( r[0] + r[1] );
  double window = r[1] - r[0];
  
  this->SetColorLevel( level );
  this->SetColorWindow( window );
}

void vtkImageView::GetColorRange( double r[2] )
{
  r[0] = this->GetColorLevel() - 0.5 * this->GetColorWindow();
  r[1] = this->GetColorLevel() + 0.5 * this->GetColorWindow();
}

//----------------------------------------------------------------------------
void vtkImageView::SetTextProperty (vtkTextProperty* textproperty)
{
  vtkSetObjectBodyMacro (TextProperty, vtkTextProperty, textproperty);
  this->CornerAnnotation->SetTextProperty (this->TextProperty);
}

//----------------------------------------------------------------------------
void vtkImageView::GetWorldCoordinatesFromImageCoordinates(int indices[3], double* position)
{
  if (!this->GetInput())
  {
    position[0] = 0; position[1] = 0; position[2] = 0;
    return;
  }
  
  // Get information
  double* spacing = this->GetInput()->GetSpacing();
  double* origin = this->GetInput()->GetOrigin();
  
  double orientedposition[4];
  for (unsigned int i=0; i<3; i++)
    orientedposition[i] = origin[i] + spacing[i]*indices[i];
  orientedposition[3] = 1;
  
  this->GetOrientationMatrix()->MultiplyPoint (orientedposition, orientedposition);
  for( unsigned int i=0; i<3; i++)
    position[i] = orientedposition[i];
}

//----------------------------------------------------------------------------
void vtkImageView::GetImageCoordinatesFromWorldCoordinates(double position[3], int* indices)
{
  if (!this->GetInput())
  {
    indices[0] = 0; indices[1] = 0; indices[2] = 0;
    return;
  }
  
  // Get information
  double unorientedposition[4] = {position[0], position[1], position[2], 1};
  double* spacing = this->GetInput()->GetSpacing();
  double* origin = this->GetInput()->GetOrigin();
  
  // apply inverted orientation matrix to the world-coordinate position
  this->InvertOrientationMatrix->MultiplyPoint (unorientedposition, unorientedposition);
  
  for (unsigned int i=0; i<3;i++)
  {
    if (fabs (spacing[i]) > 1e-15)
      indices[i] = vtkMath::Round((unorientedposition[i]-origin[i])/spacing[i]);
    else
      indices[i] = 0;
  }
}

//----------------------------------------------------------------------------
double vtkImageView::GetValueAtPosition(double worldcoordinates[3], int component )
{
  if (!this->GetInput())
    return 0.0;
  
  int indices[3];
  this->GetImageCoordinatesFromWorldCoordinates (worldcoordinates, indices);
  this->GetInputAlgorithm()->UpdateInformation();
  
  int* w_extent = this->GetInputAlgorithm()->GetOutputInformation(0)->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
  if ( (indices[0] < w_extent[0]) ||
       (indices[0] > w_extent[1]) ||
       (indices[1] < w_extent[2]) ||
       (indices[1] > w_extent[3]) ||
       (indices[2] < w_extent[4]) ||
       (indices[2] > w_extent[5]) )
    return 0;
  
  // Is the requested point in the currently loaded data extent? If not, attempt to update.
  int* extent = this->GetInput()->GetExtent();
  if ( (indices[0] < extent[0]) ||
       (indices[0] > extent[1]) ||
       (indices[1] < extent[2]) ||
       (indices[1] > extent[3]) ||
       (indices[2] < extent[4]) ||
       (indices[2] > extent[5]) ) 
  {
    int* u_extent = this->WindowLevel->GetUpdateExtent ();
    if ( (indices[0] < u_extent[0]) ||
	 (indices[0] > u_extent[1]) ||
	 (indices[1] < u_extent[2]) ||
	 (indices[1] > u_extent[3]) ||
	 (indices[2] < u_extent[4]) ||
	 (indices[2] > u_extent[5]) ) 
    {
      
      int pointExtent [6] = { indices [0], indices [0], indices [1], indices [1], indices [2], indices [2] };
      this->WindowLevel->SetUpdateExtent(pointExtent);
      this->WindowLevel->Update();
    }
    else
    {
      this->WindowLevel->Update ();
      int* new_extent = this->WindowLevel->GetOutputInformation(0)->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
      if ( (indices[0] < new_extent[0]) ||
          (indices[0] > new_extent[1]) ||
          (indices[1] < new_extent[2]) ||
          (indices[1] > new_extent[3]) ||
          (indices[2] < new_extent[4]) ||
          (indices[2] > new_extent[5]) ) 
      {
        vtkErrorMacro( "data not in slice extent after update" );
      }      
    }
  }
  else
  {  
    // Need to be sure that the input is up to date. Otherwise we may be requesting bad data.
    this->WindowLevel->Update();
  }
  
  return this->GetInput()->GetScalarComponentAsDouble (indices[0], indices[1], indices[2], component);  
}

//----------------------------------------------------------------------------
void vtkImageView::SetPosition(int a,int b) 
{
  if (this->RenderWindow)
    this->RenderWindow->SetPosition(a, b);
}

//----------------------------------------------------------------------------
int* vtkImageView::GetPosition() const
{
  if (this->RenderWindow)  
    return this->RenderWindow->GetPosition();
  return NULL;
}

//----------------------------------------------------------------------------
void vtkImageView::SetSize(int a,int b) 
{
  if (this->RenderWindow)
    this->RenderWindow->SetSize(a, b);
}

//----------------------------------------------------------------------------
int* vtkImageView::GetSize() const
{
  if (!this->RenderWindow)
    return this->RenderWindow->GetSize();
  return NULL;
}

//----------------------------------------------------------------------------
void vtkImageView::Enable (void)
{
  if (this->Interactor)
    this->Interactor->Enable();
}
//----------------------------------------------------------------------------
void vtkImageView::Disable (void)
{
  if (this->Interactor)
    this->Interactor->Disable();
}
//----------------------------------------------------------------------------
int vtkImageView::GetEnabled (void) const
{
  if (this->Interactor)
    return this->Interactor->GetEnabled();
  return -1;
}

//----------------------------------------------------------------------------
void vtkImageView::Start (void)
{
  if (this->Interactor)
    this->Interactor->Start();
}

//----------------------------------------------------------------------------
void vtkImageView::SetBackground(double rgb[3])
{
  if (this->Renderer)
    this->Renderer->SetBackground(rgb);
}
//----------------------------------------------------------------------------
double* vtkImageView::GetBackground() const
{
  if (this->Renderer)
    return this->Renderer->GetBackground();
  return NULL;
}

//----------------------------------------------------------------------------
void vtkImageView::SetZoom (double arg)
{
  if (!this->GetInput())
    return;
  
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return;  
  
  this->GetInputAlgorithm()->UpdateInformation();
  int* extent = this->GetInputAlgorithm()->GetOutputInformation(0)->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
  
  double* spacing = this->GetInput()->GetSpacing();
  double xyz[3] = {0,0,0};  
  for (unsigned int i=0; i<3; i++)
    xyz[i] = (extent [2*i +1] - extent [2*i]) * spacing[i] / 2.0;  
  double val = std::max (std::max (xyz[0], xyz[1]), xyz[2]);
  
  // Just in case no data, avoid setting scale to zero.
  val = ( val == 0. ) ? 1. : val;
	
  cam->SetParallelScale (val / arg);
  
  this->InvokeEvent (vtkImageView::ZoomChangedEvent);
  this->Modified();
}

//----------------------------------------------------------------------------
double vtkImageView::GetZoom (void)
{
  if (!this->GetInput())
    return 1.0;
  if (!this->GetInputAlgorithm() ||
      !this->GetInputAlgorithm()->GetOutputInformation(0))
    return 1.0;
  
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return 1.0;
  
  this->GetInputAlgorithm()->UpdateInformation();
  int* extent = this->GetInputAlgorithm()->GetOutputInformation(0)->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
  
  double* spacing = this->GetInput()->GetSpacing();
  double xyz[3] = {0,0,0};  
  for (unsigned int i=0; i<3; i++)
    xyz[i] = (extent [2*i +1] - extent [2*i]) * spacing[i] / 2.0;
  //xyz[i] = dimensions[i] * spacing[i] / 2.0;
  double val = std::max (std::max (xyz[0], xyz[1]), xyz[2]);
  
  return (val / cam->GetParallelScale());
}

//----------------------------------------------------------------------------
void vtkImageView::ResetCamera (void)
{
  if (this->Renderer)
  {
     if ( this->GetInput () ) 
      {
          double bounds [6];
          this->GetInputBoundsInWorldCoordinates (bounds);			  
          this->Renderer->ResetCamera (bounds);
      } else {
          // No op.
          this->Renderer->ResetCamera();
      }
     this->SetZoom (1.0);
     this->InvokeEvent (vtkImageView::CameraChangedEvent);
  }
}

//----------------------------------------------------------------------------
void vtkImageView::SetCameraPosition (double* arg)
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return;
  cam->SetPosition (arg);
  this->InvokeEvent (vtkImageView::CameraChangedEvent);
  this->Modified();
}

//----------------------------------------------------------------------------
double* vtkImageView::GetCameraPosition (void) const
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return NULL;
  return cam->GetPosition ();
}

//----------------------------------------------------------------------------
void vtkImageView::SetCameraFocalPoint (double* arg)
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return;
  cam->SetFocalPoint (arg);
  this->InvokeEvent (vtkImageView::CameraChangedEvent);
  this->Modified();
}

//----------------------------------------------------------------------------
double* vtkImageView::GetCameraFocalPoint (void) const
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return NULL;
  return cam->GetFocalPoint ();
}

//----------------------------------------------------------------------------
void vtkImageView::SetCameraViewUp (double* arg)
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return;
  cam->SetViewUp (arg);
  this->InvokeEvent (vtkImageView::CameraChangedEvent);
  this->Modified();
}

//----------------------------------------------------------------------------
double* vtkImageView::GetCameraViewUp (void) const
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return NULL;
  return cam->GetViewUp ();
}


//----------------------------------------------------------------------------
void vtkImageView::SetCameraParallelScale (double arg)
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return;
  cam->SetParallelScale (arg);
  this->InvokeEvent (vtkImageView::CameraChangedEvent);
  this->Modified();
}

//----------------------------------------------------------------------------
double vtkImageView::GetCameraParallelScale (void) const
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (!cam)
    return 1.0;
  return cam->GetParallelScale ();
}

//----------------------------------------------------------------------------
void vtkImageView::Reset (void)
{
  this->ResetCurrentPoint();
  this->ResetWindowLevel();
  // this->SetColorWindow (VTK_DOUBLE_MAX); // NT: ?? --> when i press reset I would like the windowlevels to be "reset" ?
  this->ResetCamera();	
}

//----------------------------------------------------------------------------
void vtkImageView::SetShowAnnotations (int val)
{
  this->ShowAnnotations = val;
  this->CornerAnnotation->SetVisibility (val);
}

//----------------------------------------------------------------------------
void vtkImageView::SetShowScalarBar (int val)
{
  this->ShowScalarBar = val;
  this->ScalarBar->SetVisibility (val);
}

//----------------------------------------------------------------------------
void vtkImageView::ResetWindowLevel()
{
  if (!this->GetInput())
  {
    return;
  }
  
  if( this->GetInput()->GetScalarType()==VTK_UNSIGNED_CHAR  &&
     (this->GetInput()->GetNumberOfScalarComponents()==3 || this->GetInput()->GetNumberOfScalarComponents()==4) )
  {
    return;
  }
  
  double* range = this->GetInput()->GetScalarRange();
  double window = range[1]-range[0];
  double level = 0.5*(range[1]+range[0]);  
  
  this->SetColorWindow ( window );
  this->SetColorLevel ( level );
}
//----------------------------------------------------------------------------
void vtkImageView::RemoveDataSet (vtkPointSet *arg)
{
  this->DataSetActorCollection->RemoveItem (this->FindDataSetActor (arg));
  this->DataSetCollection->RemoveItem (arg);
}

//----------------------------------------------------------------------------
vtkProp3D* vtkImageView::FindDataSetActor (vtkDataSet* arg) 
{
  int id = this->DataSetCollection->IsItemPresent (arg);
  if (id == 0)
    return NULL;
  return vtkProp3D::SafeDownCast (this->DataSetActorCollection->GetItemAsObject (id-1));
}

//----------------------------------------------------------------------------
vtkDataSet* vtkImageView::FindActorDataSet (vtkProp3D* arg) 
{
  int id = this->DataSetActorCollection->IsItemPresent (arg);
  if (id == 0)
    return NULL;
  return vtkDataSet::SafeDownCast (this->DataSetCollection->GetItemAsObject (id-1));
}



/////////////////////////////////////////////////////////////////////////////
/////////////////// NOTE ON TIME HANDLING AND ITK-BRIDGE ////////////////////
/////////////////////////////////////////////////////////////////////////////
/// Nicolas Toussaint.

/// All of this pipelines of time handling, extraction etc
/// have been implemented both for images and meshes in
/// a vtk class called vtkMetaDataSetSequence.
/// Note also that in the vtkMetaDataSetSequence implementation
/// the template switch case between scalar type has ALSO been implemented,
/// as well as all setITKInput APIs and the corr. macros...
///
/// I believe that the bridge between ITK and VTK should be implemented
/// outside of the views, from an external class, such as vtkMetaDataSet(Sequence).
/// The fact that all of these are implemented here again will make difficult the
/// maintenance of both APIs...

/// time should not be handled in this class.

void vtkImageView::GetInputBounds ( double * bounds )
{
  this->GetInputAlgorithm()->UpdateInformation();
  const int* wholeExtent = this->GetInputAlgorithm()->GetOutputInformation(0)->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
  const double * spacing = this->GetInput ()->GetSpacing ();
  const double * origin = this->GetInput ()->GetOrigin ();
	
  for ( int i(0); i < 3; ++i ) 
  {
    bounds[ 2*i     ] = spacing [ i ]*wholeExtent [ 2*i     ] + origin [ i ];
    bounds[ 2*i + 1 ] = spacing [ i ]*wholeExtent [ 2*i + 1 ] + origin [ i ];		
  }  
}

void vtkImageView::GetInputBoundsInWorldCoordinates ( double * bounds )
{
  double imageBounds [6];
  
  this->GetInputAlgorithm()->UpdateInformation();
  const int* wholeExtent = this->GetInputAlgorithm()->GetOutputInformation(0)->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
  const double * spacing = this->GetInput ()->GetSpacing ();
  const double * origin = this->GetInput ()->GetOrigin ();
  
  for ( int i(0); i < 3; ++i ) 
  {
    imageBounds[ 2*i     ] = spacing [ i ]*wholeExtent [ 2*i     ] + origin [i];
    imageBounds[ 2*i + 1 ] = spacing [ i ]*wholeExtent [ 2*i + 1 ] + origin [i];		
  }
  
  // Loop over the eight points that define the external vertices of the bounding cuboid.
  double testPoint[4];
  bounds [0] = VTK_DOUBLE_MAX;
  bounds [1] = VTK_DOUBLE_MIN;
  bounds [2] = VTK_DOUBLE_MAX;
  bounds [3] = VTK_DOUBLE_MIN;
  bounds [4] = VTK_DOUBLE_MAX;
  bounds [5] = VTK_DOUBLE_MIN;
  for ( int k(0); k<2; ++k ) {
    for ( int j(0); j<2; ++j ) {
      for ( int i(0); i<2; ++i ) {
        
        testPoint [0] = imageBounds [i  ];   // x coordinate
        testPoint [1] = imageBounds [2+j];   // y coordinate
        testPoint [2] = imageBounds [4+k];   // z coordinate
        testPoint [3] = 1.;
        
        // Transform to world coordinates
        this->GetOrientationMatrix()->MultiplyPoint (testPoint, testPoint);
        
        // Compare min / max for each coordinate.
        for ( int m(0); m < 3; ++m ) {
          
          if ( bounds [2*m] > testPoint [m] )
            bounds [2*m] = testPoint [m];
          if ( bounds [2*m + 1] < testPoint [m] )
            bounds [2*m + 1] = testPoint [m];
          
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
vtkAlgorithm* vtkImageView::GetInputAlgorithm(void)
{
  return this->WindowLevel->GetInputAlgorithm();
}


//----------------------------------------------------------------------------
void vtkImageView::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  
  os << indent << "OrientationMatrix:\n";
  this->OrientationMatrix->PrintSelf(os,indent.GetNextIndent());
  
  if ( this->LookupTable != NULL )
  {
    os << indent << "LookupTable:\n";
    this->LookupTable->PrintSelf(os,indent.GetNextIndent());
  }
  
  os << indent << "WindowLevel:\n";
  this->WindowLevel->PrintSelf(os,indent.GetNextIndent());
  
  if (this->Input)
  {
    os << indent << "Input:\n";
    this->Input->PrintSelf(os,indent.GetNextIndent());
  }
  if (this->RenderWindow)
  {
    os << indent << "RenderWindow:\n";
    this->RenderWindow->PrintSelf(os,indent.GetNextIndent());
  }
  if (this->Renderer)
  {
    os << indent << "Renderer:\n";
    this->Renderer->PrintSelf(os,indent.GetNextIndent());
  }  
  if (this->InteractorStyle)
  {
    os << indent << "InteractorStyle:\n";
    this->InteractorStyle->PrintSelf(os,indent.GetNextIndent());
  }
}
