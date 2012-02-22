/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkProlateSpheroidalTransformTensorMeshFilter.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_ProlateSpheroidalTransformTensorMeshFilter_txx_
#define _itk_ProlateSpheroidalTransformTensorMeshFilter_txx_
#include "itkProlateSpheroidalTransformTensorMeshFilter.h"

#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkProgressReporter.h"

namespace itk
{


  template <class TMesh>
  void
  ProlateSpheroidalTransformTensorMeshFilter<TMesh>
  ::PrintSelf (std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf (os,indent);
  }

  template <class TMesh>
  ProlateSpheroidalTransformTensorMeshFilter<TMesh>
  ::ProlateSpheroidalTransformTensorMeshFilter()
  {
    m_Transform = TransformType::New();
    this->ReleaseDataBeforeUpdateFlagOff();
  }
  
  template <class TMesh>
  void
  ProlateSpheroidalTransformTensorMeshFilter<TMesh>
  ::SetTransform (TransformType * transform)
  {
    m_Transform = transform;
    this->Modified();
  }
  
  template <class TMesh>
  void
  ProlateSpheroidalTransformTensorMeshFilter<TMesh>
  ::GenerateData()
  {

    
    typename MeshType::ConstPointer input  = this->GetInput();
    typename MeshType::Pointer      output = this->GetOutput();

    
    typename MeshType::PointsContainer::ConstPointer inPoints  = input->GetPoints();
    typename MeshType::PointsContainer::Pointer      outPoints = output->GetPoints();
    typename MeshType::PointDataContainer::ConstPointer inData  = input->GetPointData();
    typename MeshType::PointDataContainer::Pointer      outData = output->GetPointData();
    typename MeshType::PointsContainer::ConstIterator  it_inPoints  = inPoints->Begin();
    typename MeshType::PointsContainer::Iterator       it_outPoints = outPoints->Begin();
    typename MeshType::PointDataContainer::ConstIterator  it_inData  = inData->Begin();
    typename MeshType::PointDataContainer::Iterator       it_outData = outData->Begin();
    
    // Support for progress methods/callbacks
    ProgressReporter progress(this, 0, input->GetNumberOfPoints());

    MatrixType jacobian, transpose;
    TensorType T;
    TensorType Ts;
    PointType x;
    
    
    while( it_inPoints != inPoints->End() ) 
    {
      x = it_inPoints.Value();
      it_outPoints.Value() = m_Transform->TransformPoint (x);

      if (it_inData !=  inData->End())
      {
	// jacobian is the matrix composed ny the column vectors n, a, d
	// Hence it is the rotation matrix transfomation to apply to a vector
	// express in R0 (xyz basis) to get its R2 (prolate) coordinates.
	// But tensors are express in R1 (image basis). 
	// Direction is the transformation matrix from R1 to R0
	// Hence we have Ts = T.ApplyMatrix (Direction)
	//               Ts = Ts.ApplyMatrix (R)
	jacobian = m_Transform->GetJacobianWithRespectToCoordinates (x);
	transpose = jacobian.GetTranspose();
	
	T = it_inData.Value();
	
	if ( T.IsZero() )
	{
	  Ts = T;
	}
	else
	{
	  // Ts = T.ApplyMatrix (jacobian);
	  Ts = T.ApplyMatrix (transpose);
	}
	
	it_outData.Value() = Ts;
	++it_inData; ++it_outData;
      }

      progress.CompletedPixel();
      ++it_inPoints; ++it_outPoints;
    }
    
  }

  template <class TMesh>
  void
  ProlateSpheroidalTransformTensorMeshFilter<TMesh>
  ::GenerateOutputInformation()
  {
    // call the superclass's implementation of this method
    Superclass::GenerateOutputInformation();

    if (!m_Transform)
      itkExceptionMacro(<< "prolate spheroid not set");

    // this->GenerateData();

    
    this->CopyInputMeshToOutputMeshPoints();
    this->CopyInputMeshToOutputMeshPointData();
    this->CopyInputMeshToOutputMeshCellLinks();
    this->CopyInputMeshToOutputMeshCells();
    this->CopyInputMeshToOutputMeshCellData();
    
  }
  


} // end of namespace



#endif
