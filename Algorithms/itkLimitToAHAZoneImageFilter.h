/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkLimitToAHAZoneImageFilter.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_LimitToAHAZoneImageFilter_h_
#define _itk_LimitToAHAZoneImageFilter_h_

#include <itkImageToImageFilter.h>
#include <itkProlateSpheroidalTransform.h>
#include <itkLinearInterpolateImageFunction.h>

namespace itk
{

  class Zone
  {
  public:
    
    Zone()
    {
      for (unsigned int i=0; i<3; i++)
      {
	this->Minimum[i] = vcl_numeric_limits<double>::min();
	this->Maximum[i] = vcl_numeric_limits<double>::max();
      }
    }
    ~Zone(){};

    void SetMinimum (double* p)
    {
      for (unsigned int i=0; i<3; i++)
	this->Minimum[i] = p[i];

      while (this->Minimum[2] < 0.0)
      {
	this->Minimum[2] = 2.0 * vnl_math::pi + this->Minimum[2];
      }
      while (this->Minimum[2] >= 2.0 * vnl_math::pi)
      {
	this->Minimum[2] = this->Minimum[2] - 2.0 * vnl_math::pi;
      }
            
    }
    void SetMaximum (double* p)
    {
      for (unsigned int i=0; i<3; i++)
	this->Maximum[i] = p[i];
      
      while (this->Maximum[2] <= 0.0)
      {
	this->Maximum[2] = 2.0 * vnl_math::pi + this->Maximum[2];
      }
      while (this->Maximum[2] > 2.0 * vnl_math::pi)
      {
	this->Maximum[2] = this->Maximum[2] - 2.0 * vnl_math::pi;
      }
    }
    double* GetMinimum (void) 
    { return this->Minimum; }
    double* GetMaximum (void) 
    { return this->Maximum; }
    
    bool IsInside (double p[3])
    {
      while (p[2] < 0.0)
      {
	p[2] = 2.0 * vnl_math::pi + p[2];
      }
      while (p[2] > 2.0 * vnl_math::pi)
      {
	p[2] = p[2] - 2.0 * vnl_math::pi;
      }

      for (unsigned int i=0 ;i<3; i++)
      {
	if ( (i == 2) && ( this->Minimum[2] > this->Maximum[2] ) )
	{
	  if ( (p[i] <= this->Minimum[i]) && (p[i] > this->Maximum[i]) )
	    return false;
	  else
	    continue;
	}
	
	if (p[i] <= this->Minimum[i])
	  return false;
	if (p[i] > this->Maximum[i])
	  return false;
      }
      return true;
    }
    
    void GetCenter (double* arg)
    {
      for (unsigned int i=0; i<3; i++)
      {
	/* if (Minimum[i] > Maximum[i]) */
	/*   arg[i] = (Maximum[i] - Minimum[i]) / 2.0; */
	/* else */
	  arg[i] = (Maximum[i] + Minimum[i]) / 2.0;
	  if (Minimum[i] > Maximum[i])
	    arg[i] += vnl_math::pi;
      }
    }	
    
  private:
    
    double Minimum[3];
    double Maximum[3];

  };
  
  
  template <class TImage>
    class ITK_EXPORT LimitToAHAZoneImageFilter :
  public ImageToImageFilter<TImage, TImage>
  {

  public:

    typedef TImage ImageType;
    typedef typename ImageType::PixelType PixelType;
    typedef typename ImageType::RegionType ImageRegionType;

    typedef double ScalarType;
    
    typedef LimitToAHAZoneImageFilter Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(LimitToAHAZoneImageFilter, ImageToImageFilter);

    itkStaticConstMacro( ImageDimension, unsigned int, ImageType::ImageDimension);

    typedef Vector<ScalarType, 3>              DisplacementType;
    typedef Image<DisplacementType, 3>    DisplacementFieldType;
    typedef ProlateSpheroidalTransform<ScalarType> TransformType;
    typedef typename ImageType::PointType PointType;

    typedef LinearInterpolateImageFunction<DisplacementFieldType, ScalarType> InterpolatorType;
    typedef ContinuousIndex< ScalarType, ImageDimension > ContinuousIndexType;

    enum AHASegmentationTypeIds
    {
      AHA_13_ZONES = 0,
      AHA_17_ZONES,
      AHA_25_ZONES,
      AHA_3_ZONES,
      AHA_1_ZONE
    };

    itkGetMacro (AHASegmentationType, unsigned int);
    itkSetClampMacro (AHASegmentationType, unsigned int, AHA_13_ZONES, AHA_1_ZONE);

    itkGetMacro (CanineDivisions, unsigned int);
    void SetCanineDivisions (unsigned int arg)
    { m_CanineDivisions = arg; this->Modified(); }
    itkBooleanMacro (CanineDivisions);
    
    itkGetMacro (AHAZone, unsigned int);
    void SetAHAZone (unsigned int zone)
    { m_AHAZone = zone; this->Modified();}

    void SetDisplacementField (DisplacementFieldType::Pointer field)
    {
      m_DisplacementField  = field;
      this->Modified();
    }
    void SetInverseDisplacementField (DisplacementFieldType::Pointer field)
    {
      m_InverseDisplacementField  = field;
      this->Modified();
    }
    
    void SetTransform (TransformType::Pointer forward)
    {
      m_Transform = forward;
      this->CalculateZones();
      this->Modified();
    }

    void SetVentricleSizes(double thickness, double maxangle)
    {
      m_Thickness = thickness;
      m_MaxAngle = maxangle;
      this->Modified();
    }

    PointType GetZoneCentralPointCartesian (void);
    PointType GetZoneCentralPointProlate (void);

    unsigned int GetNumberOfAHAZones (void) const
    {
      unsigned int ret = 0;
      switch(m_AHASegmentationType)
      {
	  case AHA_13_ZONES:
	    ret = 13;
	    break;
	  case AHA_17_ZONES:
	    ret = 17;
	    break;
	  case AHA_25_ZONES:
	    ret = 25;
	    break;
	  case AHA_1_ZONE:
	    ret = 1;
	    break;
	  case AHA_3_ZONES:
	    ret = 3;
	    break;
	  default:
	    ret = 0;
	    break;
      }
      return ret;
    }

    unsigned int InWhichZoneIsPoint (PointType pt);
    unsigned int InWhichZoneIsCartesianPoint (PointType pt);
    void CalculateZones (void);
    
  protected:
    LimitToAHAZoneImageFilter();
    ~LimitToAHAZoneImageFilter(){};

    void BeforeThreadedGenerateData(void);
    void AfterThreadedGenerateData(void);
    void ThreadedGenerateData(const ImageRegionType &RegionForThread, int threadId);
    void PrintSelf(std::ostream& os, Indent indent) const
    {
      Superclass::PrintSelf(os,indent);
    }

    bool IsPointInZone (PointType point, unsigned int zone);
    
    void Calculate13Zones (void);
    void Calculate17Zones (void);
    void Calculate25Zones (void);
    void Calculate3Zones (void);
    void Calculate1Zone (void);

    typename TransformType::Pointer m_Transform;
    typename DisplacementFieldType::Pointer m_DisplacementField;
    typename DisplacementFieldType::Pointer m_InverseDisplacementField;
    typename InterpolatorType::Pointer m_Interpolator;
    typename InterpolatorType::Pointer m_InverseInterpolator;

    std::vector< Zone > m_Zones;

    double m_Thickness, m_MaxAngle;

    unsigned int m_AHAZone;
    unsigned int m_AHASegmentationType;
    unsigned int m_CanineDivisions;
    
  private:

  };


} // end of namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLimitToAHAZoneImageFilter.txx"
#endif

#endif
