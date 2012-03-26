/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkTensorImageToMeshFilter.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_TensorImageToMeshFilter_txx_
#define _itk_TensorImageToMeshFilter_txx_

#include <itkImageRegionConstIterator.h>

namespace itk
{
  
  template<class TPixel, unsigned int TDimension>
  void
  TensorImageToMeshFilter<TPixel, TDimension>
  ::SetInput(unsigned int index, const ImageType* image )
  {
    
    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(index, 
                                   const_cast< ImageType *>( image ) );
    
  }
  
  template<class TPixel, unsigned int TDimension>
  TensorImageToMeshFilter<TPixel, TDimension>
  ::TensorImageToMeshFilter()
  {
    this->SetNumberOfRequiredOutputs(1);
    this->SetNumberOfRequiredInputs(1);
  }
  
  template<class TPixel, unsigned int TDimension>
  void
  TensorImageToMeshFilter<TPixel, TDimension>
  ::GenerateData(void)
  {
    typedef ImageRegionConstIterator<ImageType> InputIteratorType;
    InputIteratorType itIn (this->GetInput(),  this->GetInput()->GetLargestPossibleRegion());
    PixelType pix;
    PointType x; x[0] = x[1] = x[2];
    unsigned long counter = 0;
    while(!itIn.IsAtEnd())
    {
      pix = itIn.Get();
      this->GetInput()->TransformIndexToPhysicalPoint(itIn.GetIndex(), x);
      
      if (pix.GetTrace() > vcl_numeric_limits<ScalarType>::epsilon())
      // if (pix.GetTrace() > 0.05)
      {
	pix = pix.ApplyMatrix(this->GetInput()->GetDirection());
	this->GetOutput()->SetPoint (counter, x);
	this->GetOutput()->SetPointData (counter, pix);
	counter++;
      }
      
      ++itIn;
    }
  }
  
  template<class TPixel, unsigned int TDimension>
  void
  TensorImageToMeshFilter<TPixel, TDimension>
  ::GenerateOutputInformation()
  {
    if (!this->GetInput())  itkExceptionMacro(<< "input image not set");
    if (!this->GetOutput()) itkExceptionMacro(<< "output mesh not set");
    
    typename ImageType::SizeType size = this->GetInput()->GetLargestPossibleRegion().GetSize();
    unsigned long maxnumberofpoints = 1;
    for (unsigned int i=0; i<TDimension; i++) maxnumberofpoints *= size[i];
    
    typename MeshType::PointsContainer::Pointer    DataPoints = MeshType::PointsContainer::New();
    typename MeshType::PointDataContainer::Pointer DataPixels = MeshType::PointDataContainer::New();  

    DataPoints->Reserve (maxnumberofpoints);
    DataPixels->Reserve (maxnumberofpoints);

    this->GetOutput()->SetPoints (DataPoints);
    this->GetOutput()->SetPointData (DataPixels);
  }
  
  
}


#endif
