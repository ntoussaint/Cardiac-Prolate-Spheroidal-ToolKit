#ifndef __itkOblateSpheroidalTransformToDomainImageFilter_h
#define __itkOblateSpheroidalTransformToDomainImageFilter_h

#include <itkObjectFactory.h>
#include <itkImageToImageFilter.h>
#include <itkOblateSpheroidalTransform.h>
#include <itkImageToVTKImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkCastImageFilter.h>

#include <vtkLandmarkSegmentationController.h>

class vtkPointSet;
class vtkContourFilter;
class vtkSmoothPolyDataFilter;

namespace itk
{

  class ITK_EXPORT OblateSpheroidalTransformToDomainImageFilter:
  public ImageToImageFilter < vtkLandmarkSegmentationController::BooleanImageType, vtkLandmarkSegmentationController::BooleanImageType>
  {
    
  public:
    typedef OblateSpheroidalTransformToDomainImageFilter          Self;
    
    typedef ImageToImageFilter<vtkLandmarkSegmentationController::BooleanImageType, vtkLandmarkSegmentationController::BooleanImageType> Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;
    
    itkNewMacro(Self);
    itkTypeMacro (OblateSpheroidalTransformToDomainImageFilter, ImageToImageFilter);

    typedef vtkLandmarkSegmentationController::ScalarType ScalarType;
    typedef vtkLandmarkSegmentationController::BooleanImageType BooleanImageType;
    typedef vtkLandmarkSegmentationController::ImageType ScalarImageType;

    typedef OblateSpheroidalTransform<ScalarType> TransformType;
    
    typedef ImageToVTKImageFilter<BooleanImageType> ConverterType;
    typedef MultiplyImageFilter<ScalarImageType,ScalarImageType> MultipyFilterType;
    typedef BooleanImageType::PointType PointType;
    typedef BooleanImageType::SpacingType VectorType;
    typedef BinaryThresholdImageFilter<ScalarImageType,ScalarImageType> ThresholdFilterType;
    typedef CastImageFilter<ScalarImageType,BooleanImageType> CastFilterType;
    
    void SetInput(const BooleanImageType *image);
    void SetCroppingPlane (PointType planeorigin, VectorType planenormal);
    
    void SetTransform (TransformType::Pointer arg);
    itkGetObjectMacro (Transform, TransformType);
    
    itkSetMacro (WallThickness, double); itkGetMacro (WallThickness, double);
    
    vtkPointSet* GetContourMesh (void)
    { return this->ContourMesh; }
    
  protected:
    OblateSpheroidalTransformToDomainImageFilter();
    ~OblateSpheroidalTransformToDomainImageFilter();
    
    void GenerateData(void);
    void GenerateOutputInformation(void);
    void GenerateDomainInformation (void);
    
    void UpdateContourMesh (void);
    
  private:
    OblateSpheroidalTransformToDomainImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    CastFilterType::Pointer      m_Caster;
    TransformType::Pointer       m_Transform;
    
    vtkPointSet* ContourMesh;
    
    PointType  m_PlaneOrigin;
    VectorType m_PlaneNormal;

    double m_WallThickness;
  };
  
  
} // end namespace itk

/* #ifndef ITK_MANUAL_INSTANTIATION */
/* #include "itkOblateSpheroidalTransformToDomainImageFilter.txx" */
/* #endif */
  
#endif
