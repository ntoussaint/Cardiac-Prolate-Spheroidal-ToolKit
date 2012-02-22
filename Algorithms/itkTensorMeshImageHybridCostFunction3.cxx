/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkTensorMeshImageHybridCostFunction3.cxx 1 2010-05-21 14:00:33Z nt08 $
Language:  C++
Author:    $Author: nt08 $
Date:      $Date: 2010-05-21 14:00:33 +0000 (Fri, 21 May 2010) $
Version:   $Revision: 1 $

Copyright (c) 2010 King's College London - Division of Imaging Sciences. All rights reserved.
See Copyright.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkTensorMeshImageHybridCostFunction3_cxx
#define _itkTensorMeshImageHybridCostFunction3_cxx

#include <itkTensorMeshImageHybridCostFunction3.h>
#include <itkTensorMeshIO.h>
#include <itkTensorImageIO.h>

namespace itk
{

TensorMeshImageHybridCostFunction3
::TensorMeshImageHybridCostFunction3()
{
  m_WarperIn       = WarperType::New();
  m_WarperOut      = WarperType::New();
  m_TransformerIn  = TransformerType::New();
  m_TransformerOut = TransformerType::New();
  m_Interpolator   = Interpolator2Type::New();  

  m_InputData                 = 0;
  m_InputReference            = 0;
  m_Transform                 = 0;
  m_ForwardDisplacementField  = 0;
  m_BackwardDisplacementField = 0;

  m_InternalData      = MeshType::New();
  m_InternalReference = MeshType::New();
  m_InternalDomain    = MeshType::New();

  m_UseProlateSpheroidalCoordinates = 0;
  m_KernelType = TensorMeshImageHybridCostFunction3::Gaussian;

  for (unsigned int i=0; i<100; i++)
  {
    m_Bounds[0][i] = vcl_numeric_limits<ParametersValueType>::min();
    m_Bounds[1][i] = vcl_numeric_limits<ParametersValueType>::max();
  }

  m_Lambda = 0.1;
  m_Beta = 0.000001;

  m_AttachTermType = MeanSquareDistance;
}

TensorMeshImageHybridCostFunction3
::~TensorMeshImageHybridCostFunction3()
{
}

  
void TensorMeshImageHybridCostFunction3
::GetOutputTerms (const ParametersType & parameters, MeasureType terms[2]) const
{
  this->GetValueInternal (parameters, terms);
}

  
TensorMeshImageHybridCostFunction3::MeasureType
TensorMeshImageHybridCostFunction3
::GetValue (const ParametersType & parameters) const
{
  
  double max = 60.00;
  
  for (unsigned int i=0; i<parameters.GetSize(); i++)
  {
    if (parameters[i] < m_Bounds[0][i])
      return max;
    if (parameters[i] > m_Bounds[1][i])
      return max;
  }
  
  double terms[2];
  this->GetValueInternal (parameters, terms);
  return terms[0] + m_Lambda * terms[1];
}
  

void TensorMeshImageHybridCostFunction3
::GetValueInternal (const ParametersType & parameters, MeasureType terms[2]) const
{
  if (!m_InternalData->GetNumberOfPoints()  || !m_InternalReference->GetNumberOfPoints())   itkExceptionMacro (<<"Missing inputs : you need to feed the Input Data and the Input Reference data through SetInputs()");
  if (m_UseProlateSpheroidalCoordinates) if ( !m_ForwardDisplacementField || !m_BackwardDisplacementField || !m_Transform ) itkExceptionMacro (<<"Unconsistent parameters : using prolate spheroidal coord. induces the need of displacement fields and transform through SetDisplacementFields() and SetTransform()");
  
  double* alpha = new double[parameters.GetSize()];
  for (unsigned int i=0; i<parameters.GetSize(); i++)
    alpha[i] = parameters[i];

  m_Interpolator->SetAlpha (alpha);
  std::cout<<"interpolating."<<std::endl;
  m_Interpolator->Update();
  MeshType::Pointer interpolatoroutput = m_Interpolator->GetOutput();
  interpolatoroutput->DisconnectPipeline();
  std::cout<<"recovering output."<<std::endl;
  MeshType::Pointer interpolatoroutput_in_cartesian = this->RecoverOutput (interpolatoroutput);
  std::cout<<"estimating attached term."<<std::endl;
  MeasureType attach = 0.0;
  switch(m_AttachTermType)
  {
      case MeanSquareDistance :
	attach = this->EstimateMeanSquaredDistanceBetweenTensorFields (interpolatoroutput, m_InternalReference);
	break;
      case MeanSquareAngle :
	attach = this->EstimateMeanSquaredAngleBetweenTensorFields (interpolatoroutput, m_InternalReference);	
	break;
      default:
	break;
  }  
  std::cout<<"estimating smoothness term."<<std::endl;
  MeasureType totalvariation = this->EstimateTensorFieldTotalVariation (interpolatoroutput_in_cartesian);
  std::cout<<"attached and gradient : "<<attach <<" -- "<<totalvariation<<std::endl;

  terms[0] = attach;
  terms[1] = totalvariation;

  delete [] alpha;
}

  
TensorMeshImageHybridCostFunction3::MeasureType
TensorMeshImageHybridCostFunction3
::EstimateTensorFieldTotalVariation( const MeshType::Pointer mesh) const
{
  ImageType::Pointer image = this->CopyMeshToImage (mesh, m_InputDomain);
  
  ImageType::SizeType radius;
  radius[0] = radius[1] = radius[2] = 1;
  
  NeighborhoodIterator<ImageType> it (radius, image, image->GetLargestPossibleRegion());
  ImageType::SpacingType spacing = image->GetSpacing();
  
  MeasureType gradientmagnitude = 0;
  unsigned int counter = 0;
  
  while (!it.IsAtEnd())
  {
    if (it.GetCenterPixel().GetTrace() > 0.0)
    {
      TensorType  L = it.GetCenterPixel(), gradient (0.0);
      
      MeasureType squaredgradient = 0.0;
  
      for(unsigned int i=0; i< ImageType::ImageDimension; i++)
      {
	TensorType Ln = it.GetNext (i);
	TensorType Lmn = it.GetPrevious (i);
	
	bool isNzero = Ln.IsZero();
	bool isMNzero = Lmn.IsZero();
	
	if( !isNzero || !isMNzero )
	{
	  // neuman conditions
	  if(isNzero && !isMNzero)
	    Ln = Lmn;
	  if(isMNzero && !isNzero) 
	    Lmn = Ln;
	}
	
	TensorType T = ( Ln - Lmn ) / static_cast<MeasureType>( 2.0 * spacing[i] );
	squaredgradient += T.GetSquaredNorm();
	// gradient += ( Ln - Lmn ) / static_cast<double>( 2 * spacing[i]*spacing[i] );
      }
      
      counter++;
      gradientmagnitude += squaredgradient + this->m_Beta*this->m_Beta;
      // gradientmagnitude += gradient.GetSquaredNorm() + this->m_Beta*this->m_Beta;
    }
    ++it;    
  }
  
  gradientmagnitude /= (double)(counter);
  
  return gradientmagnitude;
}

  
TensorMeshImageHybridCostFunction3::MeasureType
TensorMeshImageHybridCostFunction3
::EstimateMeanSquaredDistanceBetweenTensorFields( const MeshType::Pointer tensors1,
										    const MeshType::Pointer tensors2 ) const
{
  TensorType T1, T2, T;  
  MeasureType d = 0.0, trace = 0.0;
  
  PixelContainer::ConstPointer pixels1 = tensors1->GetPointData();
  PixelContainer::ConstPointer pixels2 = tensors2->GetPointData();
  PixelContainer::ConstIterator it1 = pixels1->Begin();
  PixelContainer::ConstIterator it2 = pixels2->Begin();
  
  if (!pixels1 || !pixels2)
  {
    itkWarningMacro (<<"One or both of the tensor meshes to analyze is empty of any data.");
    return 50.0;
  }
  if (tensors1->GetNumberOfPoints() != tensors1->GetNumberOfPoints())
  {
    itkWarningMacro (<<"The 2 tensor field have different number of points.");
    return 50.0;
  }

  // Since we are in the log space,
  // the log euclidean square distance between 2 tensors is :
  // d^2 = Trace ( ( (Log(S1) - log(S2) )^2 )
  
  while( (it1 != pixels1->End()) && (it2 != pixels2->End()) )
  {
    T1 = it1.Value();
    T2 = it2.Value();
    // T = (Log(S1) - log(S2) )^2
    T = (T1 - T2) * (T1 - T2);
    trace = T.GetTrace();
    d += trace;
    ++it1;
    ++it2;
  }
  
  d /= (double)(tensors1->GetNumberOfPoints());
  
  return  d;
}

  
TensorMeshImageHybridCostFunction3::MeasureType
TensorMeshImageHybridCostFunction3
::EstimateMeanSquaredAngleBetweenTensorFields( const MeshType::Pointer tensors1,
					      const MeshType::Pointer tensors2 ) const
{
  TensorType T1, T2;  
  MeasureType d = 0.0, angle = 0.0;
  TensorType::VectorType v1, v2;
  
  PixelContainer::ConstPointer pixels1 = tensors1->GetPointData();
  PixelContainer::ConstPointer pixels2 = tensors2->GetPointData();
  PixelContainer::ConstIterator it1 = pixels1->Begin();
  PixelContainer::ConstIterator it2 = pixels2->Begin();
  
  if (!pixels1 || !pixels2)
  {
    itkWarningMacro (<<"One or both of the tensor meshes to analyze is empty of any data.");
    return 10.0;
  }
  if (tensors1->GetNumberOfPoints() != tensors1->GetNumberOfPoints())
  {
    itkWarningMacro (<<"The 2 tensor field have different number of points.");
    return 50.0;
  }

  std::cout<<"estimating mean angle..."<<std::flush;
  
  while( (it1 != pixels1->End()) && (it2 != pixels2->End()) )
  {
    T1 = it1.Value();    
    T2 = it2.Value();
    v1 = T1.GetEigenvector (2);
    v2 = T2.GetEigenvector (2);

    if ( (v1 * v2) < 0.0)
      v2 = -v2;
    
    angle = std::acos(v1 * v2);
    d += angle * angle;
    ++it1;
    ++it2;
  }
  
  d /= (double)(tensors1->GetNumberOfPoints());
  std::cout<<"done."<<std::endl;
  
  return  d;
}







void TensorMeshImageHybridCostFunction3
::UpdatePipeline (void)
{
  if (!m_InputData || !m_InputReference) return;
  if (m_UseProlateSpheroidalCoordinates) if ( !m_ForwardDisplacementField || !m_BackwardDisplacementField || !m_Transform ) return;
  
  std::cout<<"Updating Pipeling... "<<std::flush;
  
  std::cout<<"1:"<<std::flush;
  m_InputDomain = this->CreateDomain (m_InputReference);
  std::cout<<"2:"<<std::flush;
  ImageToMeshType::Pointer referenceimagetomesh = ImageToMeshType::New();
  referenceimagetomesh->SetInput (m_InputReference);
  referenceimagetomesh->Update();
  std::cout<<"3:"<<std::flush;
  
  if (m_UseProlateSpheroidalCoordinates)
  {
    m_WarperIn->SetDisplacementField (m_ForwardDisplacementField);
    m_WarperIn->SetInverseDisplacementField (m_BackwardDisplacementField);
    m_WarperOut->SetDisplacementField (m_BackwardDisplacementField);
    m_WarperOut->SetInverseDisplacementField (m_ForwardDisplacementField);
    m_TransformerIn->SetTransform (m_Transform);
    TransformType::Pointer BackTransform = TransformType::New();
    m_Transform->GetInverse (BackTransform);
    m_TransformerOut->SetTransform (BackTransform);
    
    std::cout<<"4:"<<std::flush;
    m_WarperIn->SetInput (m_InputData);
    m_WarperIn->Update();
    m_TransformerIn->SetInput (m_WarperIn->GetOutput());
    m_TransformerIn->Update();
    this->Copy (m_TransformerIn->GetOutput(), m_InternalData);
    
    std::cout<<"5:"<<std::flush;
    m_WarperIn->SetInput (referenceimagetomesh->GetOutput());
    m_WarperIn->Update();
    m_TransformerIn->SetInput (m_WarperIn->GetOutput());
    m_TransformerIn->Update();
    this->Copy (m_TransformerIn->GetOutput(), m_InternalReference);    
  }
  else
  {
    this->Copy (m_InputData, m_InternalData);
    this->Copy (referenceimagetomesh->GetOutput(), m_InternalReference);
  }

  std::cout<<"6:"<<std::flush;
  this->Copy (m_InternalReference, m_InternalDomain);
  
  std::cout<<"7... "<<std::flush;
  m_Interpolator->SetInput (0, m_InternalDomain);
  m_Interpolator->SetInput (1, m_InternalData);
  m_Interpolator->SetUsePiWorkAround(m_UseProlateSpheroidalCoordinates);
  
  switch(m_KernelType)
  {
      case TensorMeshImageHybridCostFunction3::BSpline:
	m_Interpolator->SetKernel (BSplineKernelFunctionType::New());
	break;
      case TensorMeshImageHybridCostFunction3::KaiserBessel:
	m_Interpolator->SetKernel (KaiserBesselKernelType::New());
	break;
      case TensorMeshImageHybridCostFunction3::Gaussian:
      default:
	m_Interpolator->SetKernel (GaussianKernelFunctionType::New());
	break;
  }
  std::cout<<"Done"<<std::endl;

  // TensorImageIO<ScalarType,3,3>::Pointer imagewriter = TensorImageIO<ScalarType,3,3>::New();
  // TensorMeshIO<ScalarType,3,3>::Pointer  meshwriter  = TensorMeshIO<ScalarType,3,3>::New();
  // meshwriter->SetInput(m_InternalData);
  // meshwriter->SetFileName ("internaldata.vtk");
  // meshwriter->Write();
  // meshwriter->SetInput(m_InternalReference);
  // meshwriter->SetFileName ("internalreference.vtk");
  // meshwriter->Write();
  // getchar();
  
  
}

TensorMeshImageHybridCostFunction3::MeshType::Pointer
TensorMeshImageHybridCostFunction3
::RecoverOutput (MeshType::Pointer input) const
{
  MeshType::Pointer output = MeshType::New();
  
  if (m_UseProlateSpheroidalCoordinates)
  {
    TransformerType::Pointer Switcher    = TransformerType::New();
    WarperType::Pointer      Warper      = WarperType::New();
    TransformType::Pointer BackTransform = TransformType::New();
    m_Transform->GetInverse (BackTransform);
    Switcher->SetTransform (BackTransform);
    Warper->SetDisplacementField (m_BackwardDisplacementField);
    Warper->SetInverseDisplacementField (m_ForwardDisplacementField);
    
    Switcher->SetInput (input);
    Switcher->Update();
    Warper->SetInput (Switcher->GetOutput());
    Warper->Update();
    this->Copy (Warper->GetOutput(), output);
  }
  else
  {
    this->Copy (input, output);
  }

  return output;
}
  

void
TensorMeshImageHybridCostFunction3
::Copy (const MeshType::Pointer input, MeshType::Pointer output) const
{

  PointContainer::Pointer outputPoints = PointContainer::New();
  const PointContainer* inputPoints = input->GetPoints();

  if( inputPoints )
  {
    outputPoints->Reserve( inputPoints->Size() );
    
    PointContainer::ConstIterator inputItr = inputPoints->Begin();
    PointContainer::ConstIterator inputEnd = inputPoints->End();
    PointContainer::Iterator outputItr = outputPoints->Begin();
    
    while( inputItr != inputEnd )
    {
      outputItr.Value() = inputItr.Value();
      ++inputItr;
      ++outputItr;
    }
    
    output->SetPoints( outputPoints );
  }
  else
  {
    return;
  }  
  
  PixelContainer::Pointer outputPointData = PixelContainer::New();
  const PixelContainer* inputPointData = input->GetPointData();

  if( inputPointData )
  {
    outputPointData->Reserve( inputPointData->Size() );

    PixelContainer::ConstIterator inputItr = inputPointData->Begin();
    PixelContainer::ConstIterator inputEnd = inputPointData->End();
    PixelContainer::Iterator outputItr = outputPointData->Begin();
    
    while( inputItr != inputEnd )
    {
      outputItr.Value() = inputItr.Value();
      ++inputItr;
      ++outputItr;
    }
    
    output->SetPointData( outputPointData );
  }

  
}

TensorMeshImageHybridCostFunction3::ImageType::Pointer
TensorMeshImageHybridCostFunction3
::CopyMeshToImage (MeshType::Pointer mesh, DomainImageType::Pointer domain) const
{
  ImageType::Pointer image = ImageType::New();
  
  image->SetRegions (domain->GetLargestPossibleRegion());
  image->SetOrigin(domain->GetOrigin());
  image->SetSpacing(domain->GetSpacing());
  image->SetDirection(domain->GetDirection());
  image->Allocate();
  image->FillBuffer (static_cast<TensorType>(0.0));
  
  PixelContainer::ConstIterator  itIn    = mesh->GetPointData()->Begin();
  ImageRegionIterator<ImageType> itOut   (image, image->GetLargestPossibleRegion());
  ImageRegionIterator<DomainImageType> itDomain(domain, domain->GetLargestPossibleRegion());
  
  while( !itOut.IsAtEnd())
  {
    if (itDomain.Get())
    {
      itOut.Set (itIn.Value());
      ++itIn;
    }
    ++itDomain;
    ++itOut;
  }

  return image;
}

TensorMeshImageHybridCostFunction3::DomainImageType::Pointer
TensorMeshImageHybridCostFunction3
::CreateDomain (ImageType::Pointer image) const
{

  DomainImageType::Pointer domain = DomainImageType::New();
  
  domain->SetRegions (image->GetLargestPossibleRegion());
  domain->SetOrigin(image->GetOrigin());
  domain->SetSpacing(image->GetSpacing());
  domain->SetDirection(image->GetDirection());
  domain->Allocate();
  
  ImageRegionIterator<DomainImageType> itOut(domain, domain->GetLargestPossibleRegion());
  ImageRegionIterator<ImageType> itIn(image, image->GetLargestPossibleRegion());
  
  while( !itOut.IsAtEnd() )
  {
    if (!itIn.Get().IsZero())
      itOut.Set (static_cast<ScalarType>(1.0));
    else
      itOut.Set (static_cast<ScalarType>(0.0));
    ++itIn;
    ++itOut;
  }

  return domain;
  
}

} // end namespace itk
#endif
