/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkLimitToAHAZoneImageFilter.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_LimitToAHAZoneImageFilter_txx_
#define _itk_LimitToAHAZoneImageFilter_txx_

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <cmath>

namespace itk
{


  template<class TImage>
  LimitToAHAZoneImageFilter<TImage>
  ::LimitToAHAZoneImageFilter()
  {

    m_AHAZone = 0;
    
    m_Transform = NULL;
    m_DisplacementField = NULL;
    m_InverseDisplacementField = NULL;
    
    m_Interpolator = InterpolatorType::New();
    m_InverseInterpolator = InterpolatorType::New();

    m_Thickness = 15.0;
    m_MaxAngle = 93 * vnl_math::pi / 180.0;
    // m_MaxAngle = 105 * vnl_math::pi / 180.0;

    m_AHASegmentationType = AHA_17_ZONES;
    m_CanineDivisions = 0;
  }

  
  template<class TImage>
  unsigned int
  LimitToAHAZoneImageFilter<TImage>
  ::InWhichZoneIsPoint(PointType pt)
  {
    for (unsigned int zone=1; zone<=this->GetNumberOfAHAZones(); zone++)
      if (this->IsPointInZone (pt, zone))
	return zone;
    return 0;
  }
  
  template<class TImage>
  unsigned int
  LimitToAHAZoneImageFilter<TImage>
  ::InWhichZoneIsCartesianPoint(PointType pt)
  {
    
    if (!m_Transform || !m_InverseDisplacementField)
    {
      itkExceptionMacro (<<"Missing inputs for the limitation: Please provide transform and displacement field.");
    }
    
    m_InverseInterpolator->SetInputImage (m_InverseDisplacementField);    
    ContinuousIndexType index;
    bool isinside = m_InverseDisplacementField->TransformPhysicalPointToContinuousIndex (pt, index);
    
    if (isinside)
    {
      DisplacementType d = m_InverseInterpolator->EvaluateAtContinuousIndex (index);
      pt += d;
      pt = m_Transform->TransformPoint (pt);
      
      return (this->InWhichZoneIsPoint (pt));
    }
    else
      return 0;
  }
  

  template<class TImage>
  void
  LimitToAHAZoneImageFilter<TImage>
  ::BeforeThreadedGenerateData()
  {
    if (!m_Transform || !m_InverseDisplacementField)
    {
      itkExceptionMacro (<<"Missing inputs for the limitation: Please provide transform and displacement field.");
    }

    m_InverseInterpolator->SetInputImage (m_InverseDisplacementField);

    this->CalculateZones();
  }

  
  
  template<class TImage>
  void
  LimitToAHAZoneImageFilter<TImage>
  ::CalculateZones (void)
  {
    switch (m_AHASegmentationType)
    {
	case AHA_13_ZONES:
	  this->Calculate13Zones();
	  break;
	case AHA_17_ZONES:
	  this->Calculate17Zones();
	  break;
	case AHA_25_ZONES:
	  this->Calculate25Zones();
	  break;
	case AHA_3_ZONES:
	  this->Calculate3Zones();
	  break;
	case AHA_1_ZONE:
	  this->Calculate1Zone();
	  break;
	default:
	  itkExceptionMacro (<<"AHA Segmentation Type not recognized : "<<m_AHASegmentationType);
	  break;
    }
  }
  
  template<class TImage>
  void
  LimitToAHAZoneImageFilter<TImage>
  ::ThreadedGenerateData(const ImageRegionType &outputRegionForThread, int threadId)
  {
    
    typedef ImageRegionConstIterator<ImageType> InputIteratorType;
    typedef ImageRegionIterator<ImageType> OutputIteratorType;
    typedef ImageRegionIterator<DisplacementFieldType> DisplacementIteratorType;
    
    OutputIteratorType itOut(this->GetOutput(), outputRegionForThread);
    InputIteratorType itIn(this->GetInput(), outputRegionForThread);
    ContinuousIndexType index;
    PointType x;
    
    while(!itOut.IsAtEnd())
    {
      this->GetOutput()->TransformIndexToPhysicalPoint (itOut.GetIndex(), x);
      bool isinside = m_InverseDisplacementField->TransformPhysicalPointToContinuousIndex (x, index);

      if (isinside)
      {
	DisplacementType d = m_InverseInterpolator->EvaluateAtContinuousIndex (index);
	x += d;
	x = m_Transform->TransformPoint (x);

	if (this->IsPointInZone (x, m_AHAZone))
	  itOut.Set (itIn.Get());
	else
	  itOut.Set (static_cast<PixelType>(0.0));
      }
      else
	itOut.Set (static_cast<PixelType>(0.0));
      
      ++itIn;
      ++itOut;
    }
  }

  
  template<class TImage>
  void
  LimitToAHAZoneImageFilter<TImage>
  ::AfterThreadedGenerateData()
  {
    
  }
  
  template<class TImage>
  bool
  LimitToAHAZoneImageFilter<TImage>
  ::IsPointInZone (PointType point, unsigned int zone)
  {
    if (zone == 0)
      return true;
    if (zone > m_Zones.size())
      return false;
    
    unsigned int zone_id = zone - 1;
    double* p = point.GetDataPointer();
    
    return m_Zones[zone_id].IsInside (p);
  }
  
  template<class TImage>
  typename LimitToAHAZoneImageFilter<TImage>::PointType
  LimitToAHAZoneImageFilter<TImage>
  ::GetZoneCentralPointCartesian (void)
  {
    
    PointType centre = this->GetZoneCentralPointProlate();
    
    // transform centre into cartesian coordinates
    typename TransformType::Pointer inversetransform = TransformType::New();
    m_Transform->GetInverse (inversetransform);
    centre = inversetransform->TransformPoint (centre);
    
    // transform back into normal geometry
    m_Interpolator->SetInputImage (m_DisplacementField);
    ContinuousIndexType index;
    bool isinside = m_DisplacementField->TransformPhysicalPointToContinuousIndex (centre, index);
    if (!isinside )
    {
      itkExceptionMacro (<<"Point is outside the region of interest");
    }
    DisplacementType d = m_Interpolator->EvaluateAtContinuousIndex (index);
    centre += d;

    // return the centre of the zone
    return centre;
  }
  
  template<class TImage>
  typename LimitToAHAZoneImageFilter<TImage>::PointType
  LimitToAHAZoneImageFilter<TImage>
  ::GetZoneCentralPointProlate (void)
  {
    
    PointType centre; centre[0] = 0; centre[1] = 0; centre[2] = 0;
    if (m_AHAZone == 0) return centre;
    
    if (!m_Transform || !m_DisplacementField)
    {
      itkExceptionMacro (<<"Missing inputs for the centre estimation: Please provide transform and displacement fieldSS.");
    }

    // get the centre coordinates in prolate
    unsigned int zone_id = m_AHAZone - 1;
    double centred[3];
    m_Zones[zone_id].GetCenter(centred);
    for (unsigned int i=0; i<3; i++)
      centre[i] = centred[i];

    return centre;
  }
  

  template<class TImage>
  void
  LimitToAHAZoneImageFilter<TImage>
  ::Calculate17Zones()
  {
    std::cout<<"calculating for 17 zones"<<std::endl;
    
    if (!m_Transform)
    {
      itkExceptionMacro (<<"Missing inputs for calculating zones: Please provide transform and displacement field.");
    }

    double mu1 = asinh ((m_Transform->GetLambda2() - m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    double mu2 = asinh ((m_Transform->GetLambda2() + m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    
    Zone z;
    
    double p1[3] = {0,0,0};
    double p2[3] = {0,0,0};

    double spxi2 = m_MaxAngle / 4.0;
    
    m_Zones.clear();
    
    for (int i=0; i<3; i++)
    {
      int circumferential_divisions = (i < 2) ? 6 : 4;
      double shift = (i < 2) ? 0 : vnl_math::pi / 12.0;

      if (m_CanineDivisions && (i < 2))
      {
	double rv_angle = 140.0 * vnl_math::pi / 180.0;	
	double lv_division = (2.0 * vnl_math::pi - rv_angle) / 4.0;
	double rv_division = rv_angle / 2.0;
	
	p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = - lv_division;
	p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = 0.0;
	z.SetMinimum(p1); z.SetMaximum(p2);
	m_Zones.push_back (z);
	
	for (int j=0; j<2; j++)
	{
	  p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = (double)(  j) * rv_division;
	  p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = (double)(j+1) * rv_division;
	  z.SetMinimum(p1); z.SetMaximum(p2);
	  m_Zones.push_back (z);
	}
	
	for (int j=0; j<3; j++)
	{
	  p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = rv_angle + (double)(  j) * lv_division;
	  p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = rv_angle + (double)(j+1) * lv_division;
	  z.SetMinimum(p1); z.SetMaximum(p2);
	  m_Zones.push_back (z);
	}
      }
      else 
      {
	for (int j=-1; j<circumferential_divisions-1; j++)
	{
	  p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = (double)(  j) * 2.0 * vnl_math::pi / (double)(circumferential_divisions) + shift;
	  p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = (double)(j+1) * 2.0 * vnl_math::pi / (double)(circumferential_divisions) + shift;
	  z.SetMinimum(p1); z.SetMaximum(p2);
	  m_Zones.push_back (z);
	}
      }
    }
    
    p1[0] = mu1; p1[1] = 0.0 ;     p1[2] = 0.0;
    p2[0] = mu2; p2[1] = spxi2;    p2[2] = 2.0 * vnl_math::pi;
    z.SetMinimum(p1); z.SetMaximum(p2);
    m_Zones.push_back (z);
    
  }

  template<class TImage>
  void
  LimitToAHAZoneImageFilter<TImage>
  ::Calculate3Zones()
  {
    std::cout<<"calculating for 3 zones"<<std::endl;
    
    if (!m_Transform)
    {
      itkExceptionMacro (<<"Missing inputs for calculating zones: Please provide transform and displacement field.");
    }

    double mu1 = asinh ((m_Transform->GetLambda2() - m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    double mu2 = asinh ((m_Transform->GetLambda2() + m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    
    Zone z;
    
    double p1[3] = {0,0,0};
    double p2[3] = {0,0,0};

    double spxi2 = m_MaxAngle / 3.0;
    
    m_Zones.clear();
    
    for (int i=0; i<3; i++)
    {
      p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = 0;
      p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = 2.0 * vnl_math::pi;
      z.SetMinimum(p1); z.SetMaximum(p2);
      m_Zones.push_back (z);
    }
  }

  
  template<class TImage>
  void
  LimitToAHAZoneImageFilter<TImage>
  ::Calculate13Zones()
  {
    std::cout<<"calculating for 13 zones"<<std::endl;
    if (!m_Transform)
    {
      itkExceptionMacro (<<"Missing inputs for calculating zones: Please provide transform and displacement field.");
    }


    double mu1 = asinh ((m_Transform->GetLambda2() - m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    double mu2 = asinh ((m_Transform->GetLambda2() + m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    
    Zone z;
    
    double p1[3] = {0,0,0};
    double p2[3] = {0,0,0};

    double spxi2 = m_MaxAngle / 4.0;
    
    m_Zones.clear();
    
    for (int i=0; i<3; i++)
    {
      int circumferential_divisions = 4;
      double shift = 0.0;

      if (m_CanineDivisions && (i < 2))
      {
	double rv_angle = 140.0 * vnl_math::pi / 180.0;	
	double lv_division = (2.0 * vnl_math::pi - rv_angle) / 2.0;
	double rv_division = rv_angle / 2.0;
	
	p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = 0.0;
	p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = lv_division;
	z.SetMinimum(p1); z.SetMaximum(p2);
	m_Zones.push_back (z);
	
	for (int j=0; j<2; j++)
	{
	  p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = - (double)(j+1) * rv_division;
	  p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = - (double)(  j) * rv_division;
	  z.SetMinimum(p1); z.SetMaximum(p2);
	  m_Zones.push_back (z);
	}
	
	for (int j=0; j<1; j++)
	{
	  p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = - rv_angle - (double)(j+1) * lv_division;
	  p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = - rv_angle - (double)(  j) * lv_division;
	  z.SetMinimum(p1); z.SetMaximum(p2);
	  m_Zones.push_back (z);
	}
      }
      else 
      {
	for (int j=-1; j<circumferential_divisions-1; j++)
	{
	  p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = - (double)(j+1) * 2.0 * vnl_math::pi / (double)(circumferential_divisions) - shift;
	  p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = - (double)(  j) * 2.0 * vnl_math::pi / (double)(circumferential_divisions) - shift;
	  z.SetMinimum(p1); z.SetMaximum(p2);
	  m_Zones.push_back (z);
	}
      }
    }
    
    p1[0] = mu1; p1[1] = 0.0 ;     p1[2] = 0.0;
    p2[0] = mu2; p2[1] = spxi2;    p2[2] = 2.0 * vnl_math::pi;
    z.SetMinimum(p1); z.SetMaximum(p2);
    m_Zones.push_back (z);

  }

  template<class TImage>
  void
  LimitToAHAZoneImageFilter<TImage>
  ::Calculate25Zones()
  {
    std::cout<<"calculating for 25 zones"<<std::endl;
    if (!m_Transform)
    {
      itkExceptionMacro (<<"Missing inputs for calculating zones: Please provide transform and displacement field.");
    }

    double mu1 = asinh ((m_Transform->GetLambda2() - m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    double mu2 = asinh ((m_Transform->GetLambda2() + m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    
    Zone z;
    
    double p1[3] = {0,0,0};
    double p2[3] = {0,0,0};

    double spxi2 = m_MaxAngle / 4.0;
    
    m_Zones.clear();
    
    for (int i=0; i<3; i++)
    {
      int circumferential_divisions = 8;
      double shift = 0.0;

      if (m_CanineDivisions && (i < 2))
      {
	double rv_angle = 140.0 * vnl_math::pi / 180.0;
	double lv_division = (2.0 * vnl_math::pi - rv_angle) / 5.0;
	double rv_division = rv_angle / 3.0;
	
	p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = 0.0;
	p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = lv_division;
	z.SetMinimum(p1); z.SetMaximum(p2);
	m_Zones.push_back (z);
	
	for (int j=0; j<3; j++)
	{
	  p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = - (double)(j+1) * rv_division;
	  p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = - (double)(  j) * rv_division;
	  z.SetMinimum(p1); z.SetMaximum(p2);
	  m_Zones.push_back (z);
	}
	
	for (int j=0; j<4; j++)
	{
	  p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = - rv_angle - (double)(j+1) * lv_division;
	  p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = - rv_angle - (double)(  j) * lv_division;
	  z.SetMinimum(p1); z.SetMaximum(p2);
	  m_Zones.push_back (z);
	}
      }
      else 
      {
	for (int j=-1; j<circumferential_divisions-1; j++)
	{
	  p1[0] = mu1; p1[1] = m_MaxAngle - (double)(i+1)*spxi2;    p1[2] = - (double)(j+1) * 2.0 * vnl_math::pi / (double)(circumferential_divisions) - shift;
	  p2[0] = mu2; p2[1] = m_MaxAngle - (double)(  i)*spxi2;    p2[2] = - (double)(  j) * 2.0 * vnl_math::pi / (double)(circumferential_divisions) - shift;
	  z.SetMinimum(p1); z.SetMaximum(p2);
	  m_Zones.push_back (z);
	}
      }
    }
    
    p1[0] = mu1; p1[1] = 0.0 ;     p1[2] = 0.0;
    p2[0] = mu2; p2[1] = spxi2;    p2[2] = 2.0 * vnl_math::pi;
    z.SetMinimum(p1); z.SetMaximum(p2);
    m_Zones.push_back (z);
    
  }
  
  

  template<class TImage>
  void
  LimitToAHAZoneImageFilter<TImage>
  ::Calculate1Zone()
  {
    std::cout<<"calculating for 1 zones"<<std::endl;
    if (!m_Transform)
    {
      itkExceptionMacro (<<"Missing inputs for calculating zones: Please provide transform and displacement field.");
    }

    double mu1 = asinh ((m_Transform->GetLambda2() - m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    double mu2 = asinh ((m_Transform->GetLambda2() + m_Thickness/2.0) / m_Transform->GetSemiFociDistance());;
    
    Zone z;
    
    double p1[3] = {0,0,0};
    double p2[3] = {0,0,0};

    m_Zones.clear();
    
    p1[0] = mu1; p1[1] = 0;             p1[2] = 0.0;
    p2[0] = mu2; p2[1] = m_MaxAngle;    p2[2] = 2.0 * vnl_math::pi;
    z.SetMinimum(p1); z.SetMaximum(p2);
    m_Zones.push_back (z);    
  }
  
}


#endif
