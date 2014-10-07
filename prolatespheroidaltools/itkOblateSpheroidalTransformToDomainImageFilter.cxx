#ifndef _itkOblateSpheroidalTransformToDomainImageFilter_txx
#define _itkOblateSpheroidalTransformToDomainImageFilter_txx

#include <itkOblateSpheroidalTransformToDomainImageFilter.h>

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
  OblateSpheroidalTransformToDomainImageFilter
  ::OblateSpheroidalTransformToDomainImageFilter()
  {
    this->SetNumberOfRequiredOutputs(1);

    m_Transform = NULL;

    m_PlaneOrigin[0] = m_PlaneOrigin[1] = m_PlaneOrigin[2] = 0.0;
    m_PlaneNormal[0] = 1.0; m_PlaneNormal[1] = 0.0; m_PlaneNormal[2] = 0.0;    

    m_WallThickness = 13;
    
    this->ContourMesh   = vtkPolyData::New();

    m_Caster = CastFilterType::New();
  }

  // ----------------------------------------------------------------------
  // Destructor
  OblateSpheroidalTransformToDomainImageFilter
  ::~OblateSpheroidalTransformToDomainImageFilter()
  {
    this->ContourMesh->Delete();
  }
  
  // ----------------------------------------------------------------------
  // Input setting
  void OblateSpheroidalTransformToDomainImageFilter
  ::SetInput(const BooleanImageType *image)
  {
    this->Superclass::SetInput (image);
    this->GenerateDomainInformation();
  }
  
  // ----------------------------------------------------------------------
  // Input setting
  void OblateSpheroidalTransformToDomainImageFilter
  ::SetCroppingPlane (PointType planeorigin, VectorType planenormal)
  {
    m_PlaneOrigin = planeorigin;
    m_PlaneNormal = planenormal;
    this->Modified();
  }
  
  // ----------------------------------------------------------------------
  // Input setting
  void OblateSpheroidalTransformToDomainImageFilter
  ::SetTransform (TransformType::Pointer arg)
  {
    m_Transform = arg;
    
    this->GenerateDomainInformation();
  }
  
  // ----------------------------------------------------------------------
  // GenerateDomainInformation
 
  void OblateSpheroidalTransformToDomainImageFilter
  ::GenerateOutputInformation(void)
  {
    this->Superclass::GenerateOutputInformation();
    if (!this->GetInput())
      return;
    BooleanImageType::Pointer output = this->GetOutput();
    output->SetOrigin (this->GetInput()->GetOrigin());
    output->SetSpacing (this->GetInput()->GetSpacing());
    output->SetDirection (this->GetInput()->GetDirection());
    output->SetRegions (this->GetInput()->GetLargestPossibleRegion());
    output->Allocate ();
  }
  
  // ----------------------------------------------------------------------
  // GenerateDomainInformation
  void
  OblateSpheroidalTransformToDomainImageFilter
  ::GenerateDomainInformation(void)
  {
    if (!this->GetInput() || !m_Transform)
      return;

    ScalarImageType::SizeType size = this->GetInput()->GetLargestPossibleRegion().GetSize();
    
    if (size[0]*size[1]*size[2] <= 0 )
      return;
        
    double thickness = m_WallThickness;
    if (thickness > m_Transform->GetLambda1())
      thickness = m_Transform->GetLambda1() - 0.1;
    
    double mu1  = asinh ((m_Transform->GetLambda1() - thickness/2.0) / m_Transform->GetSemiFociDistance());
    double mu2  = asinh ((m_Transform->GetLambda1() + thickness/2.0) / m_Transform->GetSemiFociDistance());
    double nu1  = - vnl_math::pi / 2.0;
    double nu2  =   vnl_math::pi / 2.0;
    
    BooleanImageType::ConstPointer input = this->GetInput();
    ScalarImageType::Pointer zeta1 = ScalarImageType::New();
    zeta1->SetRegions (this->GetInput()->GetLargestPossibleRegion());
    zeta1->SetSpacing(this->GetInput()->GetSpacing());
    zeta1->SetOrigin(this->GetInput()->GetOrigin());
    zeta1->SetDirection(this->GetInput()->GetDirection());
    zeta1->Allocate();
    ScalarImageType::Pointer zeta2 = ScalarImageType::New();
    zeta2->SetRegions (this->GetInput()->GetLargestPossibleRegion());
    zeta2->SetSpacing(this->GetInput()->GetSpacing());
    zeta2->SetOrigin(this->GetInput()->GetOrigin());
    zeta2->SetDirection(this->GetInput()->GetDirection());
    zeta2->Allocate();
    itk::ImageRegionIterator<ScalarImageType> itzeta1(zeta1, zeta1->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ScalarImageType> itzeta2(zeta2, zeta2->GetLargestPossibleRegion());
    PointType in;
    PointType out;
    BooleanImageType::PointType x;
    while( !itzeta1.IsAtEnd() )
    {
      zeta1->TransformIndexToPhysicalPoint(itzeta1.GetIndex(), x);      
      for (unsigned int i=0; i<3; i++)
	in[i] = x[i];
      out = m_Transform->TransformPoint (in);
      
      itzeta1.Set (static_cast<ScalarType>(out[0]));
      itzeta2.Set (static_cast<ScalarType>(out[1]));
      ++itzeta1;
      ++itzeta2;
    }
    ThresholdFilterType::Pointer threshold1 = ThresholdFilterType::New();
    threshold1->SetInput(zeta1);
    threshold1->SetInsideValue (1.0);
    threshold1->SetOutsideValue (0.0);
    threshold1->SetUpperThreshold (mu2);
    threshold1->SetLowerThreshold (mu1);
    threshold1->Update();
    ThresholdFilterType::Pointer threshold2 = ThresholdFilterType::New();
    threshold2->SetInput(zeta2);
    threshold2->SetInsideValue (1.0);
    threshold2->SetOutsideValue (0.0);
    threshold2->SetUpperThreshold (nu2);
    threshold2->SetLowerThreshold (nu1);
    threshold2->Update();
    MultipyFilterType::Pointer multiplier = MultipyFilterType::New();
    multiplier->SetInput (0, threshold1->GetOutput());
    multiplier->SetInput (1, threshold2->GetOutput());
    multiplier->Update();

    m_Caster->SetInput (multiplier->GetOutput());
    m_Caster->Update();
  }
  
  // ----------------------------------------------------------------------
  // GenerateData
  void
  OblateSpheroidalTransformToDomainImageFilter
  ::GenerateData(void)
  {
    BooleanImageType::Pointer output = this->GetOutput();
    PointType p; VectorType v;
    itk::ImageRegionIterator<BooleanImageType> it(output, output->GetLargestPossibleRegion());
    itk::ImageRegionIterator<BooleanImageType> itdomain(m_Caster->GetOutput(), m_Caster->GetOutput()->GetLargestPossibleRegion());
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
  void
  OblateSpheroidalTransformToDomainImageFilter
  ::UpdateContourMesh (void)
  {
    ConverterType::Pointer converter = ConverterType::New();
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
      BooleanImageType::DirectionType direction = this->GetOutput()->GetDirection();
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
