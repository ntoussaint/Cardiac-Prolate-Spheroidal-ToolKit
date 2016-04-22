#ifndef __vtkOblateSpheroidalTransformController_h
#define __vtkOblateSpheroidalTransformController_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkCommand.h>

#include <itkOblateSpheroidalTransform.h>
#include <itkDoubleLandmarkSetToDomainImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkOblateSpheroidalTransformToDomainImageFilter.h>

class vtkLandmarkWidget;
class vtkCollection;
class vtkOblateSpheroidalTransformController;
class vtkPolyDataMapper;
class vtkActor;
class vtkPolyData;

//BTX
class VTK_EXPORT vtkOblateSpheroidalTransformControllerCommand : public vtkCommand
{

 public:
  static vtkOblateSpheroidalTransformControllerCommand* New()
  { return new vtkOblateSpheroidalTransformControllerCommand; }
  
  virtual void Execute ( vtkObject *caller, unsigned long, void* );
  
  void SetController (vtkOblateSpheroidalTransformController* arg);
  
 protected:
  vtkOblateSpheroidalTransformControllerCommand();
  ~vtkOblateSpheroidalTransformControllerCommand();
  
 private:
  vtkOblateSpheroidalTransformController* Controller;
  
};
//ETX

class VTK_EXPORT vtkOblateSpheroidalTransformController : public vtkPolyDataAlgorithm
{
public:
  static vtkOblateSpheroidalTransformController* New();
  vtkTypeMacro(vtkOblateSpheroidalTransformController, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  typedef itk::DoubleLandmarkSetToDomainImageFilter DomainFilterType;
  typedef DomainFilterType::BooleanImageType BooleanImageType;
  typedef DomainFilterType::ImageType ImageType;
  typedef ImageType::PixelType ScalarType;
  typedef itk::OblateSpheroidalTransform<ScalarType> TransformType;
  typedef ImageType::PointType PointType;
  typedef ImageType::SpacingType VectorType;
  typedef itk::OblateSpheroidalTransformToDomainImageFilter DomainEllipsoidFilterType;
  
  virtual void SetInput(ImageType::Pointer input);
  ImageType::Pointer GetInput();
  
  void SetTransform (TransformType::Pointer arg);  
  void SetTransform (vtkPointSet* arg);
  TransformType::Pointer GetTransform (void)
  { return m_Transform; }
  
  void SetLandmarkSet1(vtkPointSet* arg)
  { m_DomainFilter->SetLandmarkSet1 (arg); }
  void SetLandmarkSet2(vtkPointSet* arg)
  { m_DomainFilter->SetLandmarkSet2 (arg); }

  vtkGetMacro (WallThickness, double); vtkSetMacro (WallThickness, double);
  
  DomainFilterType::Pointer GetDomainFilter (void)
  { return m_DomainFilter; }
  
  void SetInteractorCollection (vtkCollection* arg);
  vtkGetObjectMacro (InteractorCollection, vtkCollection);
  
  vtkGetObjectMacro (LandmarkCollection, vtkCollection);
  vtkGetObjectMacro (TotalLandmarkCollection, vtkCollection);
  
  void SetEnabled (unsigned int arg);
  vtkGetMacro (Enabled, unsigned int);
  vtkBooleanMacro (Enabled, unsigned int);
  
  void GetLandmarkSet (vtkPolyData* arg);
  BooleanImageType::Pointer GetDomain (void)
  { return m_DomainFilter->GetOutput(); }
  BooleanImageType::Pointer GetDomainEllipsoid (void)
  {
    m_DomainEllipsoidFilter->Update();
    return m_DomainEllipsoidFilter->GetOutput();
  }

  vtkPointSet* GetDomainMesh (void)
  { return m_DomainFilter->GetContourMesh(); }
  vtkPointSet* GetDomainEllipsoidMesh (void)
  { return m_DomainEllipsoidFilter->GetContourMesh(); }

  void GuessTransform (vtkPointSet* arg);
  void VerifyOrthogonality (void);
  void RefreshConstraints (void);
  vtkLandmarkWidget* AddConstraint (double* pos, int type);
  void RemoveAllConstraints (void);
  
 protected:
  vtkOblateSpheroidalTransformController();
  ~vtkOblateSpheroidalTransformController();

  int RequestData ( vtkInformation *vtkNotUsed(request),
		    vtkInformationVector **inputVector,
		    vtkInformationVector *outputVector);
  
  void LinkInteractions (void);

  void CreateGrid (TransformType::Pointer transform,
		   double mu1, double mu2,
		   double nu1, double nu2,
		   unsigned int throughwalldivisions,
		   unsigned int longdivisions,
		   unsigned int circumdivisions,
		   vtkPolyData* grid,
		   bool stay_in_oblate = 0,
		   bool imitatecircular = 0);

  
private:
  vtkOblateSpheroidalTransformController(const vtkOblateSpheroidalTransformController&);  // Not implemented.
  void operator=(const vtkOblateSpheroidalTransformController&);  // Not implemented.
  ImageType::Pointer          m_Input;
  DomainFilterType::Pointer   m_DomainFilter;
  DomainEllipsoidFilterType::Pointer m_DomainEllipsoidFilter;
  TransformType::Pointer m_Transform;
  vtkCollection*              LandmarkCollection;
  vtkCollection*              TotalLandmarkCollection;
  vtkCollection*              InteractorCollection;
  unsigned int                Enabled;
  vtkOblateSpheroidalTransformControllerCommand* Command;

  double WallThickness;

};

#endif


