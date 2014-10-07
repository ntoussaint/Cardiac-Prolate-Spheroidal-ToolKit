#ifndef __itkProlateSpheroidalTransformToDomainImageFilter_h
#define __itkProlateSpheroidalTransformToDomainImageFilter_h

#include <itkObjectFactory.h>
#include <itkImageToImageFilter.h>
#include <itkProlateSpheroidalTransform.h>
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

  class ITK_EXPORT ProlateSpheroidalTransformToDomainImageFilter:
  public ImageToImageFilter < vtkLandmarkSegmentationController::BooleanImageType, vtkLandmarkSegmentationController::BooleanImageType>
  {
    
  public:
    typedef ProlateSpheroidalTransformToDomainImageFilter          Self;
    
    typedef ImageToImageFilter<vtkLandmarkSegmentationController::BooleanImageType, vtkLandmarkSegmentationController::BooleanImageType> Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;
    
    itkNewMacro(Self);
    itkTypeMacro (ProlateSpheroidalTransformToDomainImageFilter, ImageToImageFilter);

    typedef vtkLandmarkSegmentationController::ScalarType ScalarType;
    typedef vtkLandmarkSegmentationController::BooleanImageType BooleanImageType;
    typedef vtkLandmarkSegmentationController::ImageType ScalarImageType;

    typedef ProlateSpheroidalTransform<ScalarType> TransformType;
    
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
    itkSetMacro (WallMaxAngle, double); itkGetMacro (WallMaxAngle, double);
    
    vtkPointSet* GetContourMesh (void)
    { return this->ContourMesh; }
    
  protected:
    ProlateSpheroidalTransformToDomainImageFilter();
    ~ProlateSpheroidalTransformToDomainImageFilter();
    
    void GenerateData(void);
    void GenerateOutputInformation(void);
    void GenerateDomainInformation (void);
    
    void UpdateContourMesh (void);
    
  private:
    ProlateSpheroidalTransformToDomainImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    CastFilterType::Pointer      m_Caster;
    TransformType::Pointer       m_Transform;
    
    vtkPointSet* ContourMesh;
    
    PointType  m_PlaneOrigin;
    VectorType m_PlaneNormal;

    double m_WallMaxAngle, m_WallThickness;
  };
  
  
} // end namespace itk

/* #ifndef ITK_MANUAL_INSTANTIATION */
/* #include "itkProlateSpheroidalTransformToDomainImageFilter.txx" */
/* #endif */
  
#endif
