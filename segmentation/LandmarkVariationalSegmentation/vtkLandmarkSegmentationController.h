#ifndef __vtkLandmarkSegmentationController_h
#define __vtkLandmarkSegmentationController_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkCommand.h>

#include <itkVariationalFunctionImageToImageFilter.h>
#include <itkImageToVTKImageFilter.h>

class vtkMatrixToLinearTransform;
class vtkContourFilter;
class vtkLandmarkWidget;
class vtkCollection;
class vtkLandmarkSegmentationController;
class vtkPolyDataMapper;
class vtkActor;

//BTX
class VTK_EXPORT vtkLandmarkSegmentationControllerCommand : public vtkCommand
{

 public:
  static vtkLandmarkSegmentationControllerCommand* New()
  { return new vtkLandmarkSegmentationControllerCommand; }
  
  virtual void Execute ( vtkObject *caller, unsigned long, void* );
  
  void SetController (vtkLandmarkSegmentationController* arg);
  
 protected:
  vtkLandmarkSegmentationControllerCommand();
  ~vtkLandmarkSegmentationControllerCommand();
  
 private:
  vtkLandmarkSegmentationController* Controller;
  
};
//ETX

class VTK_EXPORT vtkLandmarkSegmentationController : public vtkPolyDataAlgorithm
{
public:
  static vtkLandmarkSegmentationController* New();
  vtkTypeMacro(vtkLandmarkSegmentationController, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  typedef double ScalarType;
  typedef itk::Image<ScalarType, 3> ImageType;
  typedef itk::Image<int, 3> BooleanImageType;
  typedef itk::VariationalFunctionImageToImageFilter<ImageType> FilterType;
  typedef itk::ImageToVTKImageFilter<ImageType> ConverterType;
  typedef FilterType::ConstraintListType ConstraintListType;
  typedef FilterType::ConstraintType ConstraintType;
  
  virtual void SetInput(ImageType::Pointer input);
  ImageType::Pointer GetInput();
  
  void SetConstraints (ConstraintListType arg);  
  void SetConstraints (vtkPointSet* arg);  
  
  void SetInteractorCollection (vtkCollection* arg);
  vtkGetObjectMacro (InteractorCollection, vtkCollection);

  vtkGetObjectMacro (LandmarkCollection, vtkCollection);
  vtkGetObjectMacro (TotalLandmarkCollection, vtkCollection);
  
  void SetEnabled (unsigned int arg);
  vtkGetMacro (Enabled, unsigned int);
  vtkBooleanMacro (Enabled, unsigned int);
  
  ImageType::Pointer GetImplicitImage (void)
  { return m_Filter->GetOutput(); }
  void GetLandmarkSet (vtkPolyData* arg);
  
  void RefreshConstraints (void);
  
  vtkLandmarkWidget* AddConstraint (double* pos, int type);
  void RemoveConstraint (vtkLandmarkWidget* arg);
  
 protected:
  vtkLandmarkSegmentationController();
  ~vtkLandmarkSegmentationController();

  int RequestData ( vtkInformation *vtkNotUsed(request),
		    vtkInformationVector **inputVector,
		    vtkInformationVector *outputVector);
  
  void LinkInteractions (void);
  
private:
  vtkLandmarkSegmentationController(const vtkLandmarkSegmentationController&);  // Not implemented.
  void operator=(const vtkLandmarkSegmentationController&);  // Not implemented.
  ImageType::Pointer          m_Input;
  FilterType::Pointer         m_Filter;
  ConverterType::Pointer      m_Converter;
  vtkMatrixToLinearTransform* Transformer;
  vtkContourFilter*           SurfaceExtractor;
  ConstraintListType          Constraints;
  vtkCollection*              LandmarkCollection;
  vtkCollection*              TotalLandmarkCollection;
  vtkCollection*              InteractorCollection;
  unsigned int                Enabled;
  vtkLandmarkSegmentationControllerCommand* Command;
  vtkPolyDataMapper*          Mapper;
  vtkActor*                   Actor;
  double                      LandmarkRadius;

};

#endif


