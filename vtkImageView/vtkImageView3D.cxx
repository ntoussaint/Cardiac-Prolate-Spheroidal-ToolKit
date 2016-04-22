#include "vtkImageView3D.h"

#ifndef VTK_MAJOR_VERSION
#  include "vtkVersion.h"
#endif

#include <vtkBoundingBox.h>
#include <vtkInteractorStyleTrackball.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>
#include <vtkImageData.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkProperty.h>
#include <vtkVolume.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageActor.h>
#include <vtkAxes.h>
#include <vtkMatrix4x4.h>
#include <vtkTubeFilter.h>
#include <vtkLookupTable.h>
#include <vtkPropAssembly.h>
#include <vtkAxesActor.h>
#include <vtkTextProperty.h>
#include <vtkCaptionActor2D.h>
#include <vtkPointData.h>
#include <vtkImageBlend.h>
#include <vtkImageReslice.h>
#include "vtkRenderWindow.h"
#include "vtkScalarsToColors.h"
#include <vtkCamera.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include <vtkGeometryFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkOrientedBoxWidget.h>
#include <vtkPlaneWidget.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkVolumeProperty.h>
#include <vtkImageMapToColors.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetCollection.h>
#include <vtkProp3DCollection.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkImageView3DCroppingBoxCallback.h>
#include <vtkObjectFactory.h>


vtkStandardNewMacro(vtkImageView3D);

//----------------------------------------------------------------------------
vtkImageView3D::vtkImageView3D()
{
  this->VolumeProperty = vtkVolumeProperty::New();
  this->VolumeActor    = vtkVolume::New();
  this->VolumeMapper   = vtkSmartVolumeMapper::New();
  this->OpacityTransferFunction = vtkPiecewiseFunction::New();
  this->ColorTransferFunction = vtkColorTransferFunction::New();
  this->Callback       = vtkImageView3DCroppingBoxCallback::New();
  this->BoxWidget      = vtkOrientedBoxWidget::New();
  this->PlaneWidget    = vtkPlaneWidget::New();
  this->Cube           = vtkAnnotatedCubeActor::New();
  this->ExtraPlaneCollection      = vtkProp3DCollection::New();
  this->ExtraPlaneInputCollection = vtkProp3DCollection::New();
  
  this->SetupVolumeRendering();
  this->SetupWidgets();  
  
  this->ShowAnnotationsOn();
  this->TextProperty->SetColor (0,0,0);
  double white[3] = {0.9, 0.9, 0.9};
  this->SetBackground (white);
  
  this->RenderingMode = PLANAR_RENDERING;
  this->CroppingMode  = CROPPING_OFF;
  
  vtkInteractorStyleSwitch* styleswitch = vtkInteractorStyleSwitch::New();
  styleswitch->SetCurrentStyleToTrackballCamera();
  this->SetInteractorStyle ( styleswitch );
  styleswitch->Delete();
}

//----------------------------------------------------------------------------
vtkImageView3D::~vtkImageView3D()
{
  this->VolumeMapper->Delete();
  this->VolumeProperty->Delete();
  this->VolumeActor->Delete();
  this->OpacityTransferFunction->Delete();
  this->ColorTransferFunction->Delete();
  this->BoxWidget->Delete();
  this->Callback->Delete();
  this->Cube->Delete();
  this->PlaneWidget->Delete();
  this->ExtraPlaneCollection->Delete();
  this->ExtraPlaneInputCollection->Delete();
}

//----------------------------------------------------------------------------
unsigned long vtkImageView3D::GetMTime()
{
  typedef unsigned long MTimeType;
  
  MTimeType mTime = Superclass::GetMTime();
  
  vtkObject * objectsToInclude[] = {
    this->VolumeMapper,
    this->VolumeProperty,
    this->VolumeActor,
    this->BoxWidget,
    this->Cube,
    this->PlaneWidget
  };
  
  const int numObjects = sizeof(objectsToInclude) / sizeof(vtkObject *);
  
  for ( int i(0); i<numObjects; ++i ) {
    vtkObject * object = objectsToInclude[i];
    if (object) {
      const MTimeType testMtime = object->GetMTime();
      if ( testMtime > mTime )
	mTime = testMtime;
    }
  }
  
  this->ExtraPlaneCollection->InitTraversal();
  vtkObject* item = this->ExtraPlaneCollection->GetNextProp3D();
  while(item)
  {
    const MTimeType testMtime = item->GetMTime();
    if ( testMtime > mTime )
      mTime = testMtime;
    item = this->ExtraPlaneCollection->GetNextProp3D();
  }
  
  return mTime;
}


//----------------------------------------------------------------------------
void vtkImageView3D::ResetCamera (void)
{
  if (this->Renderer)
  {
    this->Renderer->ResetCamera();
    this->SetZoom (1.0);
    this->InvokeEvent (vtkImageView::CameraChangedEvent);
  }
}


//----------------------------------------------------------------------------
void vtkImageView3D::SetupVolumeRendering()
{
  this->SetCroppingModeToInside();

  this->VolumeProperty->SetColor        (this->ColorTransferFunction );
  this->VolumeProperty->SetScalarOpacity(this->OpacityTransferFunction );
  this->VolumeProperty->IndependentComponentsOn();
  this->VolumeProperty->SetInterpolationTypeToLinear();
  this->VolumeProperty->ShadeOff();
  this->VolumeProperty->SetDiffuse (0.9);
  this->VolumeProperty->SetAmbient (0.15);
  this->VolumeProperty->SetSpecular (0.3);
  this->VolumeProperty->SetSpecularPower (15.0);
  
  this->VolumeActor->SetProperty ( this->VolumeProperty );
  this->VolumeActor->SetMapper   ( this->VolumeMapper );
  this->VolumeActor->PickableOff();
  this->VolumeActor->DragableOff();  
  this->VolumeActor->SetVisibility (0);

  this->Callback->SetVolumeMapper ( this->VolumeMapper );
}

//----------------------------------------------------------------------------
void vtkImageView3D::SetupWidgets()
{
  // Create an annotated cube actor (directions)
  this->Cube->SetXPlusFaceText ("L");
  this->Cube->SetXMinusFaceText ("R");
  this->Cube->SetYPlusFaceText ("P");
  this->Cube->SetYMinusFaceText ("A");
  this->Cube->SetZPlusFaceText ("S");
  this->Cube->SetZMinusFaceText ("I");
  this->Cube->SetZFaceTextRotation (90);
  this->Cube->SetFaceTextScale (0.65);
  this->Cube->GetCubeProperty()->SetColor (0.5, 1, 1);
  this->Cube->GetTextEdgesProperty()->SetLineWidth (1);
  this->Cube->GetTextEdgesProperty()->SetDiffuse (0);
  this->Cube->GetTextEdgesProperty()->SetAmbient (1);
  this->Cube->GetTextEdgesProperty()->SetColor (0.18, 0.28, 0.23);
  
  this->Cube->SetTextEdgesVisibility (1);  
  this->Cube->SetCubeVisibility(1);  
  this->Cube->SetFaceTextVisibility(1);
  
  this->Cube->GetXPlusFaceProperty()->SetColor (1, 0, 0);
  this->Cube->GetXPlusFaceProperty()->SetInterpolationToFlat();
  this->Cube->GetXMinusFaceProperty()->SetColor (1, 0, 0);
  this->Cube->GetXMinusFaceProperty()->SetInterpolationToFlat();
  this->Cube->GetYPlusFaceProperty()->SetColor (0, 1, 0);
  this->Cube->GetYPlusFaceProperty()->SetInterpolationToFlat();
  this->Cube->GetYMinusFaceProperty()->SetColor (0, 1, 0);
  this->Cube->GetYMinusFaceProperty()->SetInterpolationToFlat();
  this->Cube->GetZPlusFaceProperty()->SetColor (0, 0, 1);
  this->Cube->GetZPlusFaceProperty()->SetInterpolationToFlat();
  this->Cube->GetZMinusFaceProperty()->SetColor (0, 0, 1);
  this->Cube->GetZMinusFaceProperty()->SetInterpolationToFlat();
  
  this->BoxWidget->RotationEnabledOff();
  this->BoxWidget->SetPlaceFactor (0.5);
  this->BoxWidget->SetKeyPressActivationValue ('b');
  this->BoxWidget->AddObserver (vtkCommand::InteractionEvent, this->Callback);  
  
  this->PlaneWidget->SetKeyPressActivationValue ('p');
}

//----------------------------------------------------------------------------
void vtkImageView3D::InstallPipeline()
{
  this->Superclass::InstallPipeline();
  
  if (this->Renderer)
    this->Renderer->AddViewProp (this->VolumeActor);
}

//----------------------------------------------------------------------------
void vtkImageView3D::UnInstallPipeline()
{  
  if (this->Renderer)
    this->Renderer->RemoveViewProp (this->VolumeActor);  
  this->Superclass::UnInstallPipeline();
  this->IsInteractorInstalled = 0;
}

//----------------------------------------------------------------------------
void vtkImageView3D::InstallInteractor()
{
  if (this->Interactor && this->InteractorStyle)
    this->Interactor->SetInteractorStyle (this->InteractorStyle);
  if (this->Interactor && this->RenderWindow)
  {
    this->Interactor->SetRenderWindow(this->RenderWindow);
    this->BoxWidget->SetInteractor ( this->Interactor );
    this->PlaneWidget->SetInteractor ( this->Interactor );
  }	
  this->IsInteractorInstalled = 1;
}

//----------------------------------------------------------------------------
void vtkImageView3D::UnInstallInteractor()
{
  this->BoxWidget->SetInteractor (NULL);
  this->PlaneWidget->SetInteractor (NULL);
  
  if (this->Interactor)
  {
    this->Interactor->SetRenderWindow (NULL);
    this->Interactor->SetInteractorStyle (NULL);
  }
  this->IsInteractorInstalled = 0;    
}

//----------------------------------------------------------------------------
void vtkImageView3D::SetInput(vtkImageData* image, vtkMatrix4x4 *matrix)
{
  this->Superclass::SetInput (image, matrix);
  
  if( !image )
    return;
  
  // Get whole extent : More reliable than dimensions if not up-to-date.
  int * w_extent = image->GetExtent ();
  
  int size [3] = { w_extent [1] - w_extent[0], w_extent [3] - w_extent[2], w_extent [5] - w_extent[4] };
  
  if ( (size[0] < 2) || (size[1] < 2) || (size[2] < 2) )
  {
    vtkWarningMacro ( <<"Cannot do volume rendering for a single slice, skipping"<<endl);
    
    this->VolumeMapper->SetInputData( static_cast< vtkImageData * >( NULL ) );
    this->BoxWidget->SetInputData ( (vtkImageData*)0 );
    this->PlaneWidget->SetInputData ( (vtkImageData*)0 );
  }
  else
  {
    
  
  this->VolumeMapper->SetInputData (image);
  
  // If an image is already of type unsigned char, there is no need to
  // map it through a lookup table
  if ( image->GetScalarType() == VTK_UNSIGNED_CHAR &&
       (image->GetNumberOfScalarComponents() == 3 ||
	image->GetNumberOfScalarComponents() == 4 ) )
    this->VolumeProperty->IndependentComponentsOff();
  else
    this->VolumeProperty->IndependentComponentsOn();

  // Read bounds and use these to place widget, rather than force whole dataset to be read.
  double bounds [6];
  this->GetInputBounds (bounds);
  
  this->BoxWidget->SetInputData (image);
  this->BoxWidget->PlaceWidget (bounds);
  this->Callback->Execute (this->BoxWidget, 0, bounds);
  
  this->PlaneWidget->SetInputData (image);
  this->PlaneWidget->PlaceWidget(bounds);
  }
}

//----------------------------------------------------------------------------
void vtkImageView3D::SetOrientationMatrix (vtkMatrix4x4* matrix)
{
  this->Superclass::SetOrientationMatrix (matrix);
  this->VolumeActor->SetUserMatrix (this->OrientationMatrix);
  this->BoxWidget->SetOrientationMatrix (this->OrientationMatrix);
}

//----------------------------------------------------------------------------
void vtkImageView3D::SetLookupTable (vtkLookupTable* lookuptable)
{
  this->Superclass::SetLookupTable (lookuptable);  
  this->UpdateVolumeFunctions();
}

//----------------------------------------------------------------------------
void vtkImageView3D::SetColorWindow (double s)
{
  this->Superclass::SetColorWindow (s);
  this->UpdateVolumeFunctions();
}

//----------------------------------------------------------------------------
void vtkImageView3D::SetColorLevel (double s)
{
  this->Superclass::SetColorLevel (s);
  this->UpdateVolumeFunctions();
}

//----------------------------------------------------------------------------
void vtkImageView3D::UpdateVolumeFunctions()
{
  if (this->GetLookupTable() == NULL )
    return;
  
  vtkColorTransferFunction * color   =
  this->VolumeProperty->GetRGBTransferFunction();
  vtkPiecewiseFunction     * opacity =
  this->VolumeProperty->GetScalarOpacity();
  
  double colorValue[6]   = { 0.0, 0.0, 0.0, 0.0, 0.5, 0.0 };
  double opacityValue[4] = { 0.0, 0.0,           0.5, 0.0 };
  
  const double * range = this->LookupTable->GetRange();
  double width = range[1] - range[0];
  
  int numColors = this->GetLookupTable()->GetNumberOfTableValues();
  double factor = 1.0 / static_cast< double >( numColors - 1 );
  if ( color->GetSize() == numColors && opacity->GetSize() == numColors )
  {
    for( int i = 0; i < numColors; ++i )
    {
      double x = range[0] + factor * static_cast< double >( i ) * width;
      
      double * val = this->GetLookupTable()->GetTableValue( i );
      colorValue[0] = x;
      colorValue[1] = val[0];
      colorValue[2] = val[1];
      colorValue[3] = val[2];
      color->SetNodeValue( i, colorValue );
      
      opacityValue[0] = x;
      opacityValue[1] = val[3];
      opacity->SetNodeValue( i, opacityValue);
    }
  }
  else
  {
    color->RemoveAllPoints();
    opacity->RemoveAllPoints();
    
    // this->OpacityFunction->AddPoint (0.0,  0.0);
    for ( int i = 0; i < numColors; ++i )
    {
      double x = range[0] + factor * static_cast< double >( i ) * width;
      
      double * val = this->GetLookupTable()->GetTableValue( i );
      color->AddRGBPoint( x, val[0], val[1], val[2]);
      opacity->AddPoint( x, val[3] );
    }
  }
}

//----------------------------------------------------------------------------
void vtkImageView3D::SetRenderingMode(int arg)
{
  this->RenderingMode = arg;
  this->VolumeActor->SetVisibility (arg == vtkImageView3D::VOLUME_RENDERING);

  this->ExtraPlaneCollection->InitTraversal();
  vtkProp3D* item = this->ExtraPlaneCollection->GetNextProp3D();
  while(item)
  {
    item->SetVisibility(arg == vtkImageView3D::PLANAR_RENDERING);
    item = this->ExtraPlaneCollection->GetNextProp3D();
  }
}

//---------------------------------------------------------------------------
void vtkImageView3D::SetCroppingModeToOff( void )
{
  this->SetCroppingMode( CROPPING_OFF );
}

//---------------------------------------------------------------------------
void vtkImageView3D::SetCroppingModeToInside( void )
{
  this->SetCroppingMode( CROPPING_INSIDE );
}

//---------------------------------------------------------------------------
void vtkImageView3D::SetCroppingModeToOutside( void )
{
  this->SetCroppingMode( CROPPING_OUTSIDE );
}

//---------------------------------------------------------------------------
void vtkImageView3D::SetCroppingMode( unsigned int mode )
{
  this->CroppingMode = mode;
  
  switch ( mode )
  {
    case CROPPING_OFF:
      this->VolumeMapper->CroppingOff();
      
      break;
    case CROPPING_OUTSIDE:
      this->VolumeMapper->CroppingOn();
      this->VolumeMapper->SetCroppingRegionFlagsToSubVolume();
      
      break;
    case CROPPING_INSIDE:         // fall through to default
    default:                      // default is CROPPING_INSIDE
      this->VolumeMapper->CroppingOn();
      this->VolumeMapper->SetCroppingRegionFlags( 0x7ffdfff );
  }
}

//---------------------------------------------------------------------------
unsigned int vtkImageView3D::GetCroppingMode()
{
  return this->CroppingMode;
}

//----------------------------------------------------------------------------
vtkActor* vtkImageView3D::AddDataSet (vtkPointSet* arg, vtkProperty* prop)
{
  vtkDataSetSurfaceFilter* geometryextractor = vtkDataSetSurfaceFilter::New();
  vtkPolyDataNormals* normalextractor = vtkPolyDataNormals::New();
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  vtkActor* actor = vtkActor::New();
  
  normalextractor->SetFeatureAngle (90);
  ///\todo try to skip the normal extraction filter in order to
  // enhance the visualization speed when the data is time sequence.
  geometryextractor->SetInputData (arg);
  normalextractor->SetInputConnection (geometryextractor->GetOutputPort());
  mapper->SetInputConnection (normalextractor->GetOutputPort());
  mapper->ScalarVisibilityOff();
  
  actor->SetMapper (mapper);
  if (prop)
    actor->SetProperty (prop);
  else
  {
    vtkProperty* newprop = vtkProperty::New();
    newprop->SetColor (1,1,1);
    actor->SetProperty (newprop);
    newprop->Delete();
  }
  
  this->Renderer->AddViewProp(actor);
  
  mapper->Delete();
  normalextractor->Delete();
  geometryextractor->Delete();
  actor->Delete();

  // If this is the first widget to be added, reset camera
  if ( ! this->GetInput() ) {

      vtkBoundingBox box;
      box.AddBounds( arg->GetBounds() );

    double center[3];
    box.GetCenter(center);
    this->SetCurrentPoint(center);
    double bounds[6];
    box.GetBounds(bounds);
    this->Renderer->ResetCamera(bounds);

  }
  
  this->DataSetCollection->AddItem (arg);
  this->DataSetActorCollection->AddItem ( actor);
  
  // the actor is actually not deleted as it has
  // been referenced in the renderer, so we can
  // safely return it. well hopefully.
  return actor;
}

//----------------------------------------------------------------------------
void vtkImageView3D::RemoveDataSet(vtkPointSet* arg)
{
  vtkProp3D* actor = this->FindDataSetActor (arg);
  if (actor)
    this->Renderer->RemoveViewProp (actor);

  this->Superclass::RemoveDataSet (arg);
}

//----------------------------------------------------------------------------
class ImageActorCallback : public vtkCommand
{
public:
  static ImageActorCallback *New()
  { return new ImageActorCallback; }
  
  void Execute(vtkObject *caller, 
               unsigned long event, 
               void *vtkNotUsed(callData))
  {
    if ( this->Actor && this->Viewer)
    {
      vtkImageActor* imagecaller = vtkImageActor::SafeDownCast (caller);
      if (imagecaller && (event == vtkCommand::ModifiedEvent))
      {
	this->Actor->SetDisplayExtent (imagecaller->GetDisplayExtent());	
	this->Actor->SetInputData(imagecaller->GetInput());
	this->Actor->SetInterpolate(imagecaller->GetInterpolate());
	this->Actor->SetOpacity(imagecaller->GetOpacity());
	this->Actor->SetUserMatrix (imagecaller->GetUserMatrix());
	// this->Actor->SetVisibility(imagecaller->GetVisibility());
	this->Actor->Modified();

      }
      vtkImageData* image = vtkImageData::SafeDownCast (caller);
      if (image && (event == vtkCommand::ModifiedEvent))
      {
	this->Actor->Modified();
      }
    }
  }

  void SetActor (vtkImageActor* arg)
  {
    if (this->Actor == arg)
      return;
    if (this->Actor)
      this->Actor->UnRegister (this);
    this->Actor = arg;
    if (this->Actor)
      this->Actor->Register(this);
  }
  
  vtkImageActor* GetActor()
  {
    return this->Actor;
  }

protected:
  ImageActorCallback()
  {
    this->Actor = NULL;
  }
  ~ImageActorCallback()
  {
    if (this->Actor)
      this->Actor->Delete ();
  }

private:
  
  vtkImageActor* Actor;
  vtkImageView3D* Viewer;
};


//----------------------------------------------------------------------------
void vtkImageView3D::AddExtraPlane (vtkImageActor* input)
{

  if (!this->GetRenderer())
    return;
  
  ImageActorCallback* cbk = ImageActorCallback::New();
  vtkImageActor* actor = vtkImageActor::New();
  cbk->SetActor (actor);
  actor->SetInputData (input->GetInput());
  actor->SetDisplayExtent (input->GetDisplayExtent());
  actor->SetUserMatrix (input->GetUserMatrix());
  actor->SetInterpolate(input->GetInterpolate());
  actor->SetOpacity(input->GetOpacity());
  actor->SetVisibility (input->GetVisibility());
  
  input->AddObserver (vtkCommand::ModifiedEvent, cbk);
  // if (input->GetInput())
  //   input->GetInput()->AddObserver (vtkCommand::ModifiedEvent, cbk);
  
  this->GetRenderer()->AddViewProp (actor);
  this->ExtraPlaneCollection->AddItem (actor);
  this->ExtraPlaneInputCollection->AddItem (input);
  
  
  actor->Delete();
  cbk->Delete();

  /**
     IMPORTANT NOTE
     
     Adding a 2D actor in the 3D scene should be as simple as the next line
     instead of the code above...
     
     Unfortunately it does not seem to work properly. But this is something
     we should investigate in because it would be much simpler
  */
//  this->GetRenderer()->AddActor (input);
  
}

//----------------------------------------------------------------------------
void vtkImageView3D::RemoveExtraPlane (vtkImageActor* input)
{
  if (!this->GetRenderer())
    return;
  this->ExtraPlaneCollection->InitTraversal();
  this->ExtraPlaneInputCollection->InitTraversal();
  vtkProp3D* item = this->ExtraPlaneCollection->GetNextProp3D();
  vtkProp3D* iteminput = this->ExtraPlaneInputCollection->GetNextProp3D();
  while(item && iteminput)
  {
    if ( iteminput == input)
    {
      this->GetRenderer()->RemoveViewProp (item);
      break;
    }
    item = this->ExtraPlaneCollection->GetNextProp3D();
    iteminput = this->ExtraPlaneCollection->GetNextProp3D();
  }
}
