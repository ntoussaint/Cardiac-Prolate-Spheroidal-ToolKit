/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkSliceToVolumeMutualInformationCostFunction.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itkSliceToVolumeMutualInformationCostFunction_cxx
#define _itkSliceToVolumeMutualInformationCostFunction_cxx

#include "itkSliceToVolumeMutualInformationCostFunction.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

namespace itk
{
template <class TVolumeImage, class TSliceImage>
SliceToVolumeMutualInformationCostFunction<TVolumeImage, TSliceImage>
::SliceToVolumeMutualInformationCostFunction()
{
  m_Resampler = ResamplerType::New();
  m_Metric = MetricType::New();
  m_VolumeInterpolator = VolumeInterpolatorType::New();
  m_SliceInterpolator = SliceInterpolatorType::New();
  m_Input1 = SliceImageType::New();
  m_Input2 = SliceImageType::New();
  m_Transform = TransformType::New();
  m_InverseTransform = TransformType::New();
  m_NullTransform = SliceTransformType::New();

  m_SliceImage = 0;
  m_VolumeImage = 0;
  
  m_Resampler->SetInterpolator (m_VolumeInterpolator);
  m_Metric->SetInterpolator (m_SliceInterpolator);
  m_Metric->SetFixedImage (m_Input1);
  m_Metric->SetMovingImage (m_Input2);
  m_Metric->SetTransform (m_NullTransform);
  unsigned int numberOfBins = 50;
  unsigned int numberOfSamples = 10000;
  m_Metric->SetNumberOfHistogramBins( numberOfBins ); 
  m_Metric->SetNumberOfSpatialSamples( numberOfSamples );

  
}
  
template <class TVolumeImage, class TSliceImage>
SliceToVolumeMutualInformationCostFunction<TVolumeImage, TSliceImage>
::~SliceToVolumeMutualInformationCostFunction()
{
}

template <class TVolumeImage, class TSliceImage>
typename SliceToVolumeMutualInformationCostFunction<TVolumeImage, TSliceImage>
::MeasureType
SliceToVolumeMutualInformationCostFunction<TVolumeImage, TSliceImage>
::GetValue( const ParametersType & parameters ) const
{

  if (!m_SliceImage || !m_VolumeImage)
  {
    itkExceptionMacro ( << " cost function input not set !\n");
  }
  
  
  /**
     @param parameters are of dimension 3. They correspond to a 3D translation
  */
  m_Transform->SetParameters (parameters);
  m_Transform->GetInverse (m_InverseTransform);

  PointType origin = m_SliceImage->GetOrigin();
  PointType neworigin = m_Transform->TransformPoint (origin);
  
  m_Resampler->SetInput ( m_VolumeImage );
  m_Resampler->SetOutputOrigin( neworigin );
  m_Resampler->SetOutputSpacing( m_SliceImage->GetSpacing() );
  m_Resampler->SetOutputDirection (m_SliceImage->GetDirection());
  m_Resampler->SetSize(m_SliceImage->GetLargestPossibleRegion().GetSize());
  m_Resampler->Update();

  SliceSizeType size;
  for (unsigned int i=0; i<size.GetSizeDimension(); i++)
    size[i] = m_SliceImage->GetLargestPossibleRegion().GetSize()[i];
  SliceRegionType region;
  region.SetSize (size);
  
  m_Input1->Initialize();
  m_Input2->Initialize();
  m_Input1->SetRegions (region);
  m_Input2->SetRegions (region);
  m_Input1->Allocate();
  m_Input2->Allocate();
  m_Input1->FillBuffer (1.0);
  m_Input2->FillBuffer (2.0);
  
  m_Input1->SetPixelContainer(m_SliceImage->GetPixelContainer());
  m_Input2->SetPixelContainer(m_Resampler->GetOutput()->GetPixelContainer());
  
  m_Metric->SetFixedImage (m_Input1);
  m_Metric->SetMovingImage (m_Input2);
  m_Metric->SetFixedImageRegion (region);
  m_Metric->Initialize();
  
  ParametersType nullparams;
  nullparams.SetSize (parameters.GetSize());
  nullparams.Fill (0.0);

  MeasureType v = 0.0;
  try
  {
    v = m_Metric->GetValue(nullparams);
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    v = 0;
  }

  // MeasureType v2 = 1000.0 * std::exp (-v);
  std::cout<<parameters<<" --> "<<v<<std::endl;
  
  return v;
}

  
template <class TVolumeImage, class TSliceImage> 
void
SliceToVolumeMutualInformationCostFunction<TVolumeImage, TSliceImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk
#endif
