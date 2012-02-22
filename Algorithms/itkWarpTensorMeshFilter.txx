/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkWarpTensorMeshFilter.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef __itkWarpTensorMeshFilter_txx
#define __itkWarpTensorMeshFilter_txx
#include "itkWarpTensorMeshFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkAffineTensorTransform.h"
#include "itkProgressReporter.h"

namespace itk
{

/**
 * Default constructor.
 */
  template <class TMesh, class TDisplacementField>
  WarpTensorMeshFilter<TMesh, TDisplacementField>
::WarpTensorMeshFilter()
{
  // Setup the number of required inputs
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );

  m_EdgePaddingValue = NumericTraits<PixelType>::Zero;
  m_Jacobian = 0;
  m_ReorientationStrategy = FS;
  m_DisplacementField = 0;
  m_InverseDisplacementField = 0;
  m_InverseDisplacementInterpolator = DisplacementInterpolatorType::New();

  this->ReleaseDataBeforeUpdateFlagOff();
  
}

/**
 * Standard PrintSelf method.
 */
template <class TMesh, class TDisplacementField>
void
WarpTensorMeshFilter<TMesh, TDisplacementField>
::PrintSelf(std::ostream& os, Indent indent) const
{

  Superclass::PrintSelf(os, indent);

  os << indent << "EdgePaddingValue: "
       << static_cast<typename NumericTraits<PixelType>::PrintType>(m_EdgePaddingValue)
     << std::endl;
  
}

template <class TMesh, class TDisplacementField>
void
WarpTensorMeshFilter<TMesh, TDisplacementField>
::GenerateInputRequestedRegion()
{

  // call the superclass's implementation
  Superclass::GenerateInputRequestedRegion();

  // just propagate up the output requested region for the 
  // deformation field.
  DisplacementFieldPointer fieldPtr = this->GetDisplacementField();
  if( fieldPtr )
  {
    fieldPtr->SetRequestedRegion( fieldPtr->GetLargestPossibleRegion() );
  }
  
}


template <class TMesh, class TDisplacementField>
void
WarpTensorMeshFilter<TMesh, TDisplacementField>
::GenerateOutputInformation()
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  if (!this->GetInput())                     itkExceptionMacro(<< "input mesh not set");
  if (!this->GetOutput())                    itkExceptionMacro(<< "output mesh not set");
  if( !this->GetDisplacementField() )        itkExceptionMacro(<< "deformation field not set");
  if( !this->GetInverseDisplacementField() ) itkExceptionMacro(<< "inverse deformation field not set");

  typename MeshType::Pointer output = this->GetOutput();
  
  this->CopyInputMeshToOutputMeshPoints();
  this->CopyInputMeshToOutputMeshPointData();
  this->CopyInputMeshToOutputMeshCellLinks();
  this->CopyInputMeshToOutputMeshCells();
  this->CopyInputMeshToOutputMeshCellData();  

  // compute the Jacobian:
  typename JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
  jacobianFilter->SetInput( this->GetDisplacementField() );
  jacobianFilter->SetUseImageSpacing( true );
  try
  {
    jacobianFilter->Update();
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e;
    throw itk::ExceptionObject(__FILE__,__LINE__,"Error in WarpTensorMeshFilter::BeforeThreadedGenerateData()");
  }
  
  m_Jacobian = jacobianFilter->GetOutput();

  m_InverseDisplacementInterpolator->SetInputImage (this->GetInverseDisplacementField());

}



/**
 * Compute the output for the region specified by outputRegionForThread.
 */
template <class TMesh, class TDisplacementField>
void
WarpTensorMeshFilter<TMesh, TDisplacementField>
::GenerateData()
{
  typename MeshType::ConstPointer input  = this->GetInput();
  typename MeshType::Pointer      output = this->GetOutput();

  DisplacementFieldPointer fieldPtr = this->GetDisplacementField();
  DisplacementFieldPointer inversefieldPtr = this->GetInverseDisplacementField();
  JacobianPointer jacobianPtr = this->GetJacobian();
  RegionType region = fieldPtr->GetLargestPossibleRegion();
  // support progress methods/callbacks
  ProgressReporter progress(this, 0, input->GetNumberOfPoints());
  IndexType index;

  typename MeshType::PointsContainer::ConstPointer inPoints  = input->GetPoints();
  typename MeshType::PointsContainer::Pointer      outPoints = output->GetPoints();
  typename MeshType::PointDataContainer::ConstPointer inData  = input->GetPointData();
  typename MeshType::PointDataContainer::Pointer      outData = output->GetPointData();
  typename MeshType::PointsContainer::ConstIterator  it_inPoints  = inPoints->Begin();
  typename MeshType::PointsContainer::Iterator       it_outPoints = outPoints->Begin();
  typename MeshType::PointDataContainer::ConstIterator  it_inData  = inData->Begin();
  typename MeshType::PointDataContainer::Iterator       it_outData = outData->Begin();
  typedef AffineTensorTransform< double, VectorDimension >  TransformType;
  typename TransformType::Pointer transform = TransformType::New();
  
  while( it_inPoints != inPoints->End() ) 
  {
    PointType x = it_inPoints.Value();
    
    if (! m_Jacobian->TransformPhysicalPointToIndex (x, index))
    {
      progress.CompletedPixel();
      ++it_inPoints; ++it_outPoints;
      continue;
    }
    
    DisplacementType d = m_InverseDisplacementInterpolator->Evaluate (x);
    
    PointType P = x + d;
    
    if (it_inData !=  inData->End())
    {
      
      TensorType T = it_inData.Value();
      
      if (! m_Jacobian->TransformPhysicalPointToIndex (P, index))
      {
	++it_inPoints; ++it_outPoints;
	++it_inData; ++it_outData;
	progress.CompletedPixel();
	continue;
      }
      
      JacobianType jacobian = m_Jacobian->GetPixel (index);
      
      transform->SetMatrix( jacobian );
      switch (m_ReorientationStrategy)
      {
	  case FS:
	    T = transform->TransformTensorWithFS( T );
	    break;
	  case PPD:
	    T = transform->TransformTensorWithPPD( T );
	    break;
	  default:
	    T = transform->TransformTensorWithFS( T );
	    break;
      }

      it_outData.Value() = T;
      ++it_inData; ++it_outData;
    }

    progress.CompletedPixel();
    it_outPoints.Value() = P;
    ++it_inPoints; ++it_outPoints;
  }
}


} // end namespace itk

#endif
 
