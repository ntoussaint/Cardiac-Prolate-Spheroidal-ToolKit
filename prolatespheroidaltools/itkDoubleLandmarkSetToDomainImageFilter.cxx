#ifndef _itkDoubleLandmarkSetToDomainImageFilter_txx
#define _itkDoubleLandmarkSetToDomainImageFilter_txx

#include <itkDoubleLandmarkSetToDomainImageFilter.h>

#include <itkImageFileWriter.h>

#include <vtkPointSet.h>
#include <vtkContourFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkMatrix4x4.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkPointData.h>

namespace itk
{
  // ----------------------------------------------------------------------
  // Constructor
  DoubleLandmarkSetToDomainImageFilter
  ::DoubleLandmarkSetToDomainImageFilter()
  {
    this->SetNumberOfRequiredOutputs(1);
    
    m_Resampler     = ResampleFilterType::New();
    m_Thresholder11 = ThresholdFilterType::New();
    m_Thresholder12 = ThresholdFilterType::New();
    m_Thresholder21 = ThresholdFilterType::New();
    m_Thresholder22 = ThresholdFilterType::New();
    m_Caster1       = CastFilterType::New();
    m_Caster2       = CastFilterType::New();
    m_AndFilter     = AndImageFilterType::New();
    m_NotFilter     = NotImageFilterType::New();
    m_Multiplier    = MultipyFilterType::New();

    m_PlaneOrigin[0] = m_PlaneOrigin[1] = m_PlaneOrigin[2] = 0.0;
    m_PlaneNormal[0] = 1.0; m_PlaneNormal[1] = 0.0; m_PlaneNormal[2] = 0.0;    
    
    this->Controler1    = vtkLandmarkSegmentationController::New();
    this->Controler2    = vtkLandmarkSegmentationController::New();
    this->LandmarkSet1  = NULL;
    this->LandmarkSet2  = NULL;
    this->ContourMesh   = vtkPolyData::New();
  }

  // ----------------------------------------------------------------------
  // Destructor
  DoubleLandmarkSetToDomainImageFilter
  ::~DoubleLandmarkSetToDomainImageFilter()
  {
    this->Controler1->Delete();
    this->Controler2->Delete();
    this->ContourMesh->Delete();
  }
  
  // ----------------------------------------------------------------------
  // Input setting
  void DoubleLandmarkSetToDomainImageFilter
  ::SetInput(const InputImageType *image)
  {
    this->Superclass::SetInput (image);
    this->GenerateDomainInformation();
  }
  
  // ----------------------------------------------------------------------
  // Input setting
  void DoubleLandmarkSetToDomainImageFilter
  ::SetLandmarkSet1 (vtkPointSet* arg)
  {
    this->LandmarkSet1 = arg;
    this->GenerateDomainInformation();
  }
  
  // ----------------------------------------------------------------------
  // Input setting
  void DoubleLandmarkSetToDomainImageFilter
  ::SetLandmarkSet2 (vtkPointSet* arg)
  {
    this->LandmarkSet2 = arg;
    this->GenerateDomainInformation();
  }
  
  // ----------------------------------------------------------------------
  // Input setting
  void DoubleLandmarkSetToDomainImageFilter
  ::SetCroppingPlane (PointType planeorigin, VectorType planenormal)
  {
    m_PlaneOrigin = planeorigin;
    m_PlaneNormal = planenormal;
    this->Modified();
  }
  
  // ----------------------------------------------------------------------
  // GenerateDomainInformation
  void
  DoubleLandmarkSetToDomainImageFilter
  ::GenerateOutputInformation(void)
  {
    this->Superclass::GenerateOutputInformation();
    if (!this->GetInput() || !this->LandmarkSet1 || !this->LandmarkSet2)
      return;
    
    m_AndFilter->GetOutput()->Update();
    BooleanImageType::Pointer output = this->GetOutput();
    output->SetOrigin (m_AndFilter->GetOutput()->GetOrigin());
    output->SetSpacing (m_AndFilter->GetOutput()->GetSpacing());
    output->SetDirection (m_AndFilter->GetOutput()->GetDirection());
    output->SetRegions (m_AndFilter->GetOutput()->GetLargestPossibleRegion());
    output->Allocate ();
  }
  
  // ----------------------------------------------------------------------
  // GenerateDomainInformation
  void
  DoubleLandmarkSetToDomainImageFilter
  ::GenerateDomainInformation(void)
  {
    if (!this->GetInput() || !this->LandmarkSet1 || !this->LandmarkSet2)
      return;

    ImageType::ConstPointer input = this->GetInput();
    ImageType::SpacingType spacing = input->GetSpacing();
    double minspacing = std::min (std::min (spacing[0], spacing[1]), spacing[2]);
    ImageType::SpacingType newspacing; newspacing[0] = newspacing[1] = newspacing[2] = minspacing;
    ImageType::SizeType size = input->GetLargestPossibleRegion().GetSize();
    ImageType::SizeType newsize; 
    for (unsigned int i=0; i<3; i++) newsize[i] = (unsigned int)((double)(size[i]) * spacing[i] / newspacing[i]);

    m_Resampler->SetInput (input);
    m_Resampler->SetOutputParametersFromImage (input);
    m_Resampler->SetOutputSpacing (newspacing);
    m_Resampler->SetSize (newsize);  
    m_Resampler->Update();
    
    this->Controler1->SetInput (m_Resampler->GetOutput());
    this->Controler1->SetConstraints (this->LandmarkSet1);
    this->Controler2->SetInput (m_Resampler->GetOutput());
    this->Controler2->SetConstraints (this->LandmarkSet2);
    this->Controler1->Update();
    this->Controler2->Update();
    ImageType::Pointer f1 = this->Controler1->GetImplicitImage();
    ImageType::Pointer f2 = this->Controler2->GetImplicitImage();
    m_Thresholder11->SetInput (f1);
    m_Thresholder12->SetInput (f2);
    m_Thresholder11->ThresholdBelow (0.0);
    m_Thresholder12->ThresholdBelow (0.0);
    m_Thresholder11->SetOutsideValue (0.0);
    m_Thresholder12->SetOutsideValue (0.0);
    m_Thresholder21->SetInput (m_Thresholder11->GetOutput());
    m_Thresholder22->SetInput (m_Thresholder12->GetOutput());
    m_Thresholder21->ThresholdAbove (0.000001);
    m_Thresholder22->ThresholdAbove (0.000001);
    m_Thresholder21->SetOutsideValue (1.0);
    m_Thresholder22->SetOutsideValue (1.0);
    m_Caster1->SetInput (m_Thresholder21->GetOutput());
    m_Caster2->SetInput (m_Thresholder22->GetOutput());
    m_NotFilter->SetInput (m_Caster1->GetOutput());
    m_AndFilter->SetInput1 (m_NotFilter->GetOutput());
    m_AndFilter->SetInput2 (m_Caster2->GetOutput());
    m_AndFilter->Update();
  }
  
  // ----------------------------------------------------------------------
  // GenerateData
  void DoubleLandmarkSetToDomainImageFilter
  ::GenerateInputRequestedRegion(void)
  {
    this->Superclass::GenerateInputRequestedRegion();    
    if ( !this->GetInput() )
      return;
    InputImagePointer  inputPtr  = 
      const_cast< ImageType *>( this->GetInput() );
    ImageType::RegionType inputRegion;
    inputRegion = inputPtr->GetLargestPossibleRegion();
    inputPtr->SetRequestedRegion(inputRegion);    
    return;
  }
  
  // ----------------------------------------------------------------------
  // GenerateData
  void DoubleLandmarkSetToDomainImageFilter
  ::GenerateData(void)
  {
    BooleanImageType::Pointer output = this->GetOutput();
    PointType p; VectorType v;
    itk::ImageRegionIterator<BooleanImageType> it(output, output->GetLargestPossibleRegion());
    itk::ImageRegionIterator<BooleanImageType> itdomain(m_AndFilter->GetOutput(), m_AndFilter->GetOutput()->GetLargestPossibleRegion());
    while( !it.IsAtEnd() )
    {
      output->TransformIndexToPhysicalPoint (it.GetIndex(), p);
      v = p - m_PlaneOrigin;
      bool is_below_plane = (v * m_PlaneNormal) < 0;
      if (is_below_plane && ( itdomain.Get() > 0 ) )
    	it.Set (1);
      else
    	it.Set (0);
      ++it; ++itdomain;
    }    
    this->UpdateContourMesh();
  }


  // ----------------------------------------------------------------------
  // GenerateData
  void DoubleLandmarkSetToDomainImageFilter
  ::UpdateContourMesh (void)
  {
    typedef itk::ImageToVTKImageFilter<vtkLandmarkSegmentationController::BooleanImageType> BooleanConverterType;
    BooleanConverterType::Pointer converter = BooleanConverterType::New();
    converter->SetInput (this->GetOutput());
    converter->Update();
    vtkContourFilter* filter = vtkContourFilter::New();
    filter->SetInputData (converter->GetOutput());
    filter->SetValue (0, 0.5);
    filter->ComputeScalarsOff();
    filter->ComputeNormalsOff();
    filter->Update();
    
    if (filter->GetOutput()->GetNumberOfPoints())
    {
      vtkMatrix4x4* matrix = vtkMatrix4x4::New();
      matrix->Identity();
      ImageType::DirectionType direction = this->GetOutput()->GetDirection();
      PointType origin = this->GetOutput()->GetOrigin();
      for (unsigned int i=0; i<3; i++)
	for (unsigned int j=0; j<3; j++)
	  matrix->SetElement (i,j,direction[i][j]);
      PointType correctedorigin;
      matrix->MultiplyPoint (origin.GetDataPointer(), correctedorigin.GetDataPointer());
      for (int i=0; i<3; i++)
	matrix->SetElement (i, 3, origin[i]-correctedorigin[i]);

      
      vtkSmoothPolyDataFilter* smoother = vtkSmoothPolyDataFilter::New();
      smoother->SetInputData (filter->GetOutput());
      smoother->SetNumberOfIterations (30);
      smoother->SetRelaxationFactor (0.15);
      smoother->Update();
      vtkMatrixToLinearTransform* transform = vtkMatrixToLinearTransform::New();
      transform->SetInput (matrix);
      vtkPoints* newpoints = vtkPoints::New();
      transform->TransformPoints (smoother->GetOutput()->GetPoints(),newpoints);
      this->ContourMesh->DeepCopy (smoother->GetOutput());
      
      this->ContourMesh->SetPoints (newpoints);
      transform->Delete();
      matrix->Delete();
      smoother->Delete();
      newpoints->Delete();
    }
    filter->Delete();
  }
  
} // end namespace itk

#endif
