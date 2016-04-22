#ifndef _vtkImageView3D_h_
#define _vtkImageView3D_h_

#include "vtkImageView.h"

#include <vtkOrientedBoxWidget.h>
#include <vtkAnnotatedCubeActor.h>
#include <vtkPlaneWidget.h>
#include <vtkVolumeProperty.h>

class vtkVolumeMapper;
class vtkVolume;
class vtkPiecewiseFunction;
class vtkColorTransferFunction;
class vtkVolumeProperty;
class vtkImageActor;
class vtkAxes;
class vtkViewImage2D;
class vtkScalarsToColors;
class vtkColorTransferFunction;
class vtkImage3DDisplay;
class vtkProp3DCollection;
class vtkVolume;
class vtkImageView3DCroppingBoxCallback;

/**
   \class vtkImageView3D vtkImageView3D.h "vtkImageView3D.h"
   \brief 3D image 3D rendering viewer
   \author Pierre Fillard & Marc Traina & Nicolas Toussaint
   
   This class allows to view 3D images. Images have to be
   vtkImageData.
   volume rendering and mulptiplane reconstructions are provided
   remote plan can also be used, so can be an orientation cube, ...
*/

class VTK_EXPORT vtkImageView3D : public vtkImageView
{  
public:

  static vtkImageView3D* New();
  vtkTypeMacro(vtkImageView3D, vtkImageView);
  
  // Override vtkObject - return the maximum mtime of this and any objects owned by this.
  unsigned long GetMTime();

  // Description:
  // Rendeing Modes available.
  // PLANAR_RENDERING will render every vtkImageActor instance added with Add2DPhantom()
  // whereas VOLUME_RENDERING will render the volume added with SetInput().
  //BTX
  enum RenderingModeIds
  {
    VOLUME_RENDERING = 0,
    PLANAR_RENDERING
  };
  //ETX
  
  vtkGetObjectMacro (VolumeActor, vtkVolume);
  vtkGetObjectMacro (OpacityTransferFunction, vtkPiecewiseFunction);
  vtkGetObjectMacro (ColorTransferFunction, vtkColorTransferFunction);
  vtkGetObjectMacro (VolumeProperty, vtkVolumeProperty);
  vtkGetObjectMacro (VolumeMapper, vtkVolumeMapper);
  
  vtkGetObjectMacro (PlaneWidget, vtkPlaneWidget);
  vtkGetObjectMacro (BoxWidget, vtkOrientedBoxWidget);
  vtkGetObjectMacro (ExtraPlaneCollection, vtkProp3DCollection);  

  /** Set the box widget visibility */
  void SetShowBoxWidget (int a)
  {
    if (this->Interactor)
      this->BoxWidget->SetEnabled (a);
  }
  bool GetShowBoxWidget (void)
  {
    return this->BoxWidget->GetEnabled();
  }  
  vtkBooleanMacro (ShowBoxWidget, int);
  
  /** Set the plane widget on */
  void SetShowPlaneWidget (int a)
  {
    if (this->Interactor)
      this->PlaneWidget->SetEnabled (a);
  }
  bool GetShowPlaneWidget (void)
  {
    return this->PlaneWidget->GetEnabled();
  }
  vtkBooleanMacro (ShowPlaneWidget, int);

  /** Set the cube widget on */
  void SetShowCube (int a)
  {
    if (this->Interactor)
      this->Cube->SetVisibility (a);
  }
  bool GetShowCube (void)
  {
    return this->Cube->GetVisibility();
  }
  vtkBooleanMacro (ShowCube, int);
  
  void SetShade (int a)
  {
    this->VolumeProperty->SetShade (a);
  }
  bool GetShade (void)
  {
    return this->VolumeProperty->GetShade();
  }
  vtkBooleanMacro (Shade, int);

  /** Set the rendering mode to volume rendering (VR). */
  virtual void SetRenderingModeToVR (void)
  {this->SetRenderingMode (VOLUME_RENDERING); }
  /** Set the rendering mode to planar views. */
  virtual void SetRenderingModeToPlanar (void)
  { this->SetRenderingMode (PLANAR_RENDERING); }
  /** Set the rendering mode. */
  virtual void SetRenderingMode (int mode);
  /** Get the current rendering mode. */
  vtkGetMacro (RenderingMode, int);

  // Description:
  // Cropping Modes available.
  // CROPPING_INSIDE will crop inside the image
  // whereas CROPPING_OUTSIDE will crop outside.
  enum CroppingModeIds
  {
    CROPPING_OFF     = 0,
    CROPPING_INSIDE  = 1,
    CROPPING_OUTSIDE = 2
  };

  virtual void ResetCamera (void);

  /** Set the cropping mode */
  virtual void SetCroppingModeToOff(void);
  virtual void SetCroppingModeToInside(void);
  virtual void SetCroppingModeToOutside(void);
  virtual void SetCroppingMode(unsigned int);
  // vtkGetMacro (CroppingMode, int);
  virtual unsigned int GetCroppingMode ();
  
  virtual void SetInput (vtkImageData* input, vtkMatrix4x4 *matrix = 0);
  virtual void SetOrientationMatrix (vtkMatrix4x4* matrix);
  // Description:
  // Set window and level for mapping pixels to colors.
  virtual void SetColorWindow(double s);
  virtual void SetColorLevel(double s);
  /** Set a user-defined lookup table */
  virtual void SetLookupTable (vtkLookupTable* lookuptable);

  virtual void InstallInteractor();
  virtual void UnInstallInteractor();

  virtual vtkActor* AddDataSet (vtkPointSet* arg, vtkProperty* prop = NULL);
  virtual void RemoveDataSet (vtkPointSet* arg);

  /**
     Add an extra plane to the 3D view. the argument is an image actor
     that supposingly follows a vtkImageView2D instance. The actor will
     be displayed in the 3D scene and will be fully synchronized with
     the actor it came from.
  */
  virtual void AddExtraPlane (vtkImageActor* input); 
  virtual void RemoveExtraPlane (vtkImageActor* input); 
    
protected: 

  vtkImageView3D();
  ~vtkImageView3D();

  // Description:
  virtual void InstallPipeline();
  virtual void UnInstallPipeline();
  
  virtual void SetupVolumeRendering();
  virtual void SetupWidgets();
  virtual void UpdateVolumeFunctions();

  // volume property
  vtkVolumeProperty* VolumeProperty;
  // volume actor
  vtkVolume* VolumeActor;
  // volume mapper
  vtkVolumeMapper* VolumeMapper;
  // Opacity Function
  vtkPiecewiseFunction* OpacityTransferFunction;
  // Color Function
  vtkColorTransferFunction* ColorTransferFunction;
  
  // image 3D cropping box callback
  vtkImageView3DCroppingBoxCallback* Callback; 
  // box widget
  vtkOrientedBoxWidget* BoxWidget;
  // vtkPlane widget
  vtkPlaneWidget* PlaneWidget;
  // annotated cube actor
  vtkAnnotatedCubeActor* Cube;
  
  /**
     The ExtraPlaneCollection is a collection gathering the ImageActor
     instances that are currently displayed in addition to common ones
     (ActorX, ActorY, ActorZ).
     ExtraPlaneInputCollection is "read-only" collection to be able to
     know which inputs the ExtraPlaneCollection instances are actually
     following.
  */
  vtkProp3DCollection* ExtraPlaneCollection;
  /**
     The ExtraPlaneCollection is a collection gathering the ImageActor
     instances that are currently displayed in addition to common ones
     (ActorX, ActorY, ActorZ).
     ExtraPlaneInputCollection is "read-only" collection to be able to
     know which inputs the ExtraPlaneCollection instances are actually
     following.
  */
  vtkProp3DCollection* ExtraPlaneInputCollection;
  
  unsigned int RenderingMode;
  int          CroppingMode;
  
private:
  vtkImageView3D(const vtkImageView3D&);  // Not implemented.
  void operator=(const vtkImageView3D&);    // Not implemented.
  
};
  

#endif
