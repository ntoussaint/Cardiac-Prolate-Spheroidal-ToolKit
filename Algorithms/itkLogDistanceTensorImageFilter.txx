/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkLogDistanceTensorImageFilter.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_LogDistanceTensorImageFilter_txx_
#define _itk_LogDistanceTensorImageFilter_txx_
#include "itkLogDistanceTensorImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkProgressReporter.h"
#include <itkContinuousIndex.h>
namespace itk
{


  template<class TInputImage, class TOutputImage>
  void
  LogDistanceTensorImageFilter<TInputImage,TOutputImage>
  ::PrintSelf (std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf (os,indent);
  }



  template<class TInputImage, class TOutputImage>
  LogDistanceTensorImageFilter<TInputImage,TOutputImage>
  ::LogDistanceTensorImageFilter()
  {
    
    // this filter requires two input images
    this->SetNumberOfRequiredInputs( 2 );
    m_DistanceImage = DistanceImageType::New();
  }
  

  template<class TInputImage, class TOutputImage>
  void
  LogDistanceTensorImageFilter<TInputImage,TOutputImage>
  ::SetInput2( const TInputImage * image )
  {
    this->SetNthInput(1, const_cast<TInputImage *>( image ) );
  }


  template<class TInputImage, class TOutputImage>
  const typename LogDistanceTensorImageFilter<TInputImage,TOutputImage>
  ::InputImageType *
  LogDistanceTensorImageFilter<TInputImage,TOutputImage>
  ::GetInput2()
  {
    return static_cast< const TInputImage * >
      (this->ProcessObject::GetInput(1));
  }

  template<class TInputImage, class TOutputImage>
  void
  LogDistanceTensorImageFilter<TInputImage,TOutputImage>
  ::GenerateInputRequestedRegion() throw (InvalidRequestedRegionError)
  {
    Superclass::GenerateInputRequestedRegion();

    // this filter requires:
    // - the largeset possible region of the first image
    // - the corresponding region of the second image
    if ( this->GetInput1() )
    {
      InputImagePointer image1 =
	const_cast< InputImageType * >( this->GetInput1() );
       image1->SetRequestedRegionToLargestPossibleRegion();

      if ( this->GetInput2() )
      {
	InputImagePointer image2 =
	  const_cast< InputImageType * >( this->GetInput2() );
       image2->SetRequestedRegionToLargestPossibleRegion();
// 	image2->SetRequestedRegion( 
// 				   this->GetInput1()->GetRequestedRegion() );
      }

      DistanceImagePointerType output = const_cast< DistanceImageType * >( this->GetDistanceImage() );
      output->SetRegions (image1->GetLargestPossibleRegion());
      output->SetDirection (image1->GetDirection());
      output->SetOrigin(image1->GetOrigin());
      output->SetSpacing(image1->GetSpacing());      
      output->Allocate();
    }

    
  }
  
  template<class TInputImage, class TOutputImage>
  void
  LogDistanceTensorImageFilter<TInputImage,TOutputImage>
  ::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId)
  {

    
    // Support for progress methods/callbacks
    ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
    
    // Create an iterator that will walk the output region for this thread.
    typedef ImageRegionConstIterator<InputImageType>    InputIteratorType;
    typedef ImageRegionIterator<DistanceImageType>    OutputIteratorType;

    InputIteratorType  itIn1(this->GetInput1(), outputRegionForThread);
    InputIteratorType  itIn2(this->GetInput2(), outputRegionForThread);
    OutputIteratorType itOut(const_cast<DistanceImageType*>(this->GetDistanceImage()), outputRegionForThread);
    PointType x;
    InputPixelType T1, T2, T;
    VectorType v1, v2;
    
    DistancePixelType d;
    ContinuousIndex<RealType, ImageDimension> continuousindex;
    
    while( !itOut.IsAtEnd() )
    {
      T1 = itIn1.Get();
      T2 = itIn2.Get();

      v1 = T1.GetEigenvector (2);
      v2 = T2.GetEigenvector (2);
      v1.Normalize();
      v2.Normalize();
      
      d = static_cast<DistancePixelType>(0.0);

      if ((T1.GetTrace() > 0.05) && (T2.GetTrace() > 0.05) )
      {
	
	T = (T1 - T2);
	
	//       d = static_cast<DistancePixelType>(T.GetNorm());
	//       d *= d;
	d = std::acos (std::abs (v1*v2)) * 180.0 / vnl_math::pi;
      }
      
      itOut.Set (d);
      ++itOut;
      ++itIn1;
      ++itIn2;
      progress.CompletedPixel();
    }
  }




} // end of namespace



#endif
