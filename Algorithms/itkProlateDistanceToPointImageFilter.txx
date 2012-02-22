/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkProlateDistanceToPointImageFilter.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_ProlateDistanceToPointImageFilter_txx_
#define _itk_ProlateDistanceToPointImageFilter_txx_

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace itk
{


  template<class TImage>
  ProlateDistanceToPointImageFilter<TImage>
  ::ProlateDistanceToPointImageFilter()
  {

    m_Transform = NULL;
    m_InverseDisplacementField = NULL;
    
    m_BandwidthMatrix.SetIdentity();
    
    m_InverseInterpolator = InterpolatorType::New();
  }
    

  template<class TImage>
  void
  ProlateDistanceToPointImageFilter<TImage>
  ::BeforeThreadedGenerateData()
  {
    if (!m_Transform || !m_InverseDisplacementField)
      itkExceptionMacro (<<"Missing inputs for the limitation: Please provide transform and displacement field.");

    m_InverseInterpolator->SetInputImage (m_InverseDisplacementField);

  }
  
  template<class TImage>
  void
  ProlateDistanceToPointImageFilter<TImage>
  ::ThreadedGenerateData(const ImageRegionType &outputRegionForThread, int threadId)
  {
  
    typedef ImageRegionConstIterator<ImageType> InputIteratorType;
    typedef ImageRegionIterator<ImageType> OutputIteratorType;
    typedef ImageRegionIterator<DisplacementFieldType> DisplacementIteratorType;

    OutputIteratorType itOut(this->GetOutput(), outputRegionForThread);
    ContinuousIndexType index;
    PointType x, pointinprolate;
    DisplacementType d, dX;
    bool isinside;
    double deltaX;
    
    isinside = m_InverseDisplacementField->TransformPhysicalPointToContinuousIndex (m_Point, index);
    if (!isinside)
      itkExceptionMacro (<<"Point asked is not inside the heart.");
    d = m_InverseInterpolator->EvaluateAtContinuousIndex (index);
    pointinprolate = m_Point + d;
    pointinprolate = m_Transform->TransformPoint (pointinprolate);
    
    while(!itOut.IsAtEnd())
    {
      this->GetOutput()->TransformIndexToPhysicalPoint (itOut.GetIndex(), x);
      isinside = m_InverseDisplacementField->TransformPhysicalPointToContinuousIndex (x, index);

      if (isinside)
      {
	d = m_InverseInterpolator->EvaluateAtContinuousIndex (index);
	x += d;
	x = m_Transform->TransformPoint (x);

	dX = x - pointinprolate;
	
	while (dX[2] >= vnl_math::pi)
	  dX[2] -= 2 * vnl_math::pi;
	while (dX[2] < - vnl_math::pi)
	  dX[2] += 2 * vnl_math::pi;
	while (dX[1] >= vnl_math::pi)
	  dX[1] -= 2 * vnl_math::pi;
	while (dX[1] < - vnl_math::pi)
	  dX[1] += 2 * vnl_math::pi;

	deltaX = std::sqrt ( dX * ( m_SqInverseBandwidthMatrix * dX ) );
	//G = std::sqrt(2*vnl_math::pi) * m_Kernel->Evaluate (deltaX);
	
	itOut.Set (static_cast<PixelType>(deltaX));
      }
      else
	itOut.Set (static_cast<PixelType>(0.0));
      
      ++itOut;
    }
  }

  
  template<class TImage>
  void
  ProlateDistanceToPointImageFilter<TImage>
  ::AfterThreadedGenerateData()
  {
    
  }

}


#endif
