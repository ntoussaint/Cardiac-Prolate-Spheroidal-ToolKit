/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkTensorMeshToImageFilter.txx 1 2010-05-21 14:00:33Z nt08 $
  Language:  C++
  Date:      $Date: 2010-12-10 20:55:58 +0000 (Fri, 10 Dec 2010) $
  Version:   $Revision: 122 $
  Author:    $Author: nt08 $
  
  Copyright (c) 2010 King's College London - Division of Imaging Sciences. All rights reserved.
  See Copyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itk_TensorMeshToImageFilter_txx_
#define _itk_TensorMeshToImageFilter_txx_

#include <itkImageRegionIterator.h>

namespace itk
{
  
  template<class TPixel, unsigned int TDimension>
  void
  TensorMeshToImageFilter<TPixel, TDimension>
  ::SetInput(unsigned int index, const MeshType* mesh )
  {
    
    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(index, 
                                   const_cast< MeshType *>( mesh ) );
    
  }
  
  template<class TPixel, unsigned int TDimension>
  TensorMeshToImageFilter<TPixel, TDimension>
  ::TensorMeshToImageFilter()
  {
    this->SetNumberOfRequiredOutputs(1);
    this->SetNumberOfRequiredInputs(1);
  }
  
  template<class TPixel, unsigned int TDimension>
  void
  TensorMeshToImageFilter<TPixel, TDimension>
  ::GenerateData(void)
  {
    typedef ImageRegionIterator<ImageType> OutputIteratorType;
    typedef typename MeshType::PointType   PointType;
    
    OutputIteratorType            itOut (this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());
    typename ImageType::IndexType idOut;
    
    PixelType pix (static_cast<ScalarType>(0.0));
    PointType point; point[0] = point[1] = point[2] = 0.0;
    
    for (unsigned int i=0; i<this->GetInput()->GetNumberOfPoints(); i++)
    {
      this->GetInput()->GetPoint (i, &point);
      this->GetInput()->GetPointData (i, &pix);
      
      bool isinside = this->GetOutput()->TransformPhysicalPointToIndex (point, idOut);
      if (isinside)
      {
	itOut.SetIndex (idOut);
	pix.SetVnlMatrix(pix.ApplyMatrix (this->GetOutput()->GetDirection().GetTranspose()));
	itOut.Set (pix);
      }
    }
  }

  template<class TPixel, unsigned int TDimension>
  void
  TensorMeshToImageFilter<TPixel, TDimension>
  ::GenerateOutputInformation()
  {
    typename DomainImageType::ConstPointer domain = this->GetDomain();
    typename ImageType::Pointer            output = this->GetOutput();
    
    if( !domain)  itkExceptionMacro(<<"Missing Input domain");
    if( !output ) itkExceptionMacro(<<"Missing Output Image");
    
    output->SetRegions (domain->GetLargestPossibleRegion());
    output->SetOrigin(domain->GetOrigin());
    output->SetSpacing(domain->GetSpacing());
    output->SetDirection(domain->GetDirection());
    output->Allocate();
    
    PixelType defaultvalue( static_cast<ScalarType>(0.0) );
    output->FillBuffer (defaultvalue);
  }
  
}


#endif
