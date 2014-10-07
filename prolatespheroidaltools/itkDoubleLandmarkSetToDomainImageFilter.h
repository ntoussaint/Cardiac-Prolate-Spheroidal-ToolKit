#ifndef __itkDoubleLandmarkSetToDomainImageFilter_h
#define __itkDoubleLandmarkSetToDomainImageFilter_h

#include <itkObjectFactory.h>
#include <itkImageToImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkAndImageFilter.h>
#include <itkNotImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkMultiplyImageFilter.h>

#include <vtkLandmarkSegmentationController.h>

class vtkPointSet;
class vtkContourFilter;
class vtkSmoothPolyDataFilter;

namespace itk
{

  class ITK_EXPORT DoubleLandmarkSetToDomainImageFilter:
  public ImageToImageFilter < vtkLandmarkSegmentationController::ImageType, vtkLandmarkSegmentationController::BooleanImageType>
  {
    
  public:
    typedef DoubleLandmarkSetToDomainImageFilter          Self;
    typedef ImageToImageFilter<vtkLandmarkSegmentationController::ImageType, vtkLandmarkSegmentationController::BooleanImageType> Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;
    
    typedef vtkLandmarkSegmentationController::ImageType ImageType;
    typedef vtkLandmarkSegmentationController::BooleanImageType BooleanImageType;
    
    typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleFilterType;
    typedef itk::ThresholdImageFilter<ImageType> ThresholdFilterType;
    typedef itk::CastImageFilter<ImageType,BooleanImageType> CastFilterType;
    typedef itk::AndImageFilter<BooleanImageType> AndImageFilterType;
    typedef itk::NotImageFilter<BooleanImageType,BooleanImageType> NotImageFilterType;
    typedef itk::MultiplyImageFilter<BooleanImageType,BooleanImageType> MultipyFilterType;
    typedef itk::ImageToVTKImageFilter<BooleanImageType> ConverterType;
    typedef ImageType::PointType PointType;
    typedef ImageType::SpacingType VectorType;
    
    
    itkNewMacro(Self);
    itkTypeMacro (DoubleLandmarkSetToDomainImageFilter, ImageToImageFilter);
    
    itkStaticConstMacro(ImageDimension, unsigned int,
			BooleanImageType::ImageDimension);    

    void SetInput(const InputImageType *image);
    void SetLandmarkSet1(vtkPointSet* arg);
    void SetLandmarkSet2(vtkPointSet* arg);
    void SetCroppingPlane (PointType planeorigin, VectorType planenormal);
    
    vtkPointSet* GetContourMesh (void)
    { return this->ContourMesh; }

  protected:
    DoubleLandmarkSetToDomainImageFilter();
    ~DoubleLandmarkSetToDomainImageFilter();
    
    void GenerateData(void);
    void GenerateInputRequestedRegion( void );
    void GenerateDomainInformation(void);
    void GenerateOutputInformation(void);
    void UpdateContourMesh (void);
    
  private:
    DoubleLandmarkSetToDomainImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    vtkPointSet* LandmarkSet1;
    vtkPointSet* LandmarkSet2;
    vtkLandmarkSegmentationController* Controler1;
    vtkLandmarkSegmentationController* Controler2;
    vtkPointSet* ContourMesh;
    
    ResampleFilterType::Pointer  m_Resampler;
    ThresholdFilterType::Pointer m_Thresholder11;
    ThresholdFilterType::Pointer m_Thresholder12;
    ThresholdFilterType::Pointer m_Thresholder21;
    ThresholdFilterType::Pointer m_Thresholder22;
    CastFilterType::Pointer      m_Caster1;
    CastFilterType::Pointer      m_Caster2;
    AndImageFilterType::Pointer  m_AndFilter;
    NotImageFilterType::Pointer  m_NotFilter;
    MultipyFilterType::Pointer   m_Multiplier;
    PointType  m_PlaneOrigin;
    VectorType m_PlaneNormal;
  };
  
  
} // end namespace itk

/* #ifndef ITK_MANUAL_INSTANTIATION */
/* #include "itkDoubleLandmarkSetToDomainImageFilter.txx" */
/* #endif */
  
#endif
