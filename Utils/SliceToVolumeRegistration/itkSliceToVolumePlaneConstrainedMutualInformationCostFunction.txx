/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkSliceToVolumePlaneConstrainedMutualInformationCostFunction.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itkSliceToVolumePlaneConstrainedMutualInformationCostFunction_cxx
#define _itkSliceToVolumePlaneConstrainedMutualInformationCostFunction_cxx

#include "itkSliceToVolumePlaneConstrainedMutualInformationCostFunction.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

namespace itk
{
template <class TVolumeImage, class TSliceImage>
SliceToVolumePlaneConstrainedMutualInformationCostFunction<TVolumeImage, TSliceImage>
::SliceToVolumePlaneConstrainedMutualInformationCostFunction()
{
  m_DirectionOfConstrain.SetIdentity();
}
  
template <class TVolumeImage, class TSliceImage>
SliceToVolumePlaneConstrainedMutualInformationCostFunction<TVolumeImage, TSliceImage>
::~SliceToVolumePlaneConstrainedMutualInformationCostFunction()
{
}

template <class TVolumeImage, class TSliceImage>
typename SliceToVolumePlaneConstrainedMutualInformationCostFunction<TVolumeImage, TSliceImage>
::MeasureType
SliceToVolumePlaneConstrainedMutualInformationCostFunction<TVolumeImage, TSliceImage>
::GetValue( const ParametersType & parameters ) const
{

  /**
     @param parameters are of dimension 2. They correspond to a translation
     in the plane of m_SliceImage.

     We have to transpose this in-plane translation in a 3D translation manner,
     to feed the transform m_Transform.

     translation = parameters[0] * direction[0] + parameters[1] * direction[1]
  */

  VectorType translation;
  translation[0] = parameters[0];
  translation[1] = parameters[1];
  translation[2] = 0.0;
  DirectionType constrain = m_DirectionOfConstrain;
  
  translation = constrain * translation;

  ParametersType realparameters;
  realparameters.SetSize (3);
  for (unsigned int i=0; i<3; i++)
    realparameters[i] = translation[i];
  
  this->m_Transform->SetParameters (realparameters);
  this->m_Transform->GetInverse (this->m_InverseTransform);

  PointType origin = this->m_SliceImage->GetOrigin();
  PointType neworigin = this->m_Transform->TransformPoint (origin);
  
  this->m_Resampler->SetInput ( this->m_VolumeImage );
  this->m_Resampler->SetOutputOrigin( neworigin );
  this->m_Resampler->SetOutputSpacing( this->m_SliceImage->GetSpacing() );
  this->m_Resampler->SetOutputDirection (this->m_SliceImage->GetDirection());
  this->m_Resampler->SetSize(this->m_SliceImage->GetLargestPossibleRegion().GetSize());
  this->m_Resampler->Update();

  SliceSizeType size;
  for (unsigned int i=0; i<size.GetSizeDimension(); i++)
    size[i] = this->m_SliceImage->GetLargestPossibleRegion().GetSize()[i];
  SliceRegionType region;
  region.SetSize (size);
  
  this->m_Input1->Initialize();
  this->m_Input2->Initialize();
  this->m_Input1->SetRegions (region);
  this->m_Input2->SetRegions (region);
  this->m_Input1->Allocate();
  this->m_Input2->Allocate();
  this->m_Input1->FillBuffer (1.0);
  this->m_Input2->FillBuffer (2.0);
  
  this->m_Input1->SetPixelContainer(this->m_SliceImage->GetPixelContainer());
  this->m_Input2->SetPixelContainer(this->m_Resampler->GetOutput()->GetPixelContainer());
  
  this->m_Metric->SetFixedImage (this->m_Input1);
  this->m_Metric->SetMovingImage (this->m_Input2);
  this->m_Metric->SetFixedImageRegion (region);
  this->m_Metric->Initialize();
  
  ParametersType nullparams;
  nullparams.SetSize (parameters.GetSize());
  nullparams.Fill (0.0);

  MeasureType v = 0.0;
  try
  {
    v = this->m_Metric->GetValue(nullparams);
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
SliceToVolumePlaneConstrainedMutualInformationCostFunction<TVolumeImage, TSliceImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace itk
#endif
