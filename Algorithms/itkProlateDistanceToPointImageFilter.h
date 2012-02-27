/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkProlateDistanceToPointImageFilter.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_ProlateDistanceToPointImageFilter_h_
#define _itk_ProlateDistanceToPointImageFilter_h_

#include <itkImageToImageFilter.h>
#include <itkProlateSpheroidalTransform.h>
#include <itkLinearInterpolateImageFunction.h>

namespace itk
{

  
  template <class TImage>
    class ITK_EXPORT ProlateDistanceToPointImageFilter :
  public ImageToImageFilter<TImage, TImage>
  {

  public:

    
    typedef ProlateDistanceToPointImageFilter Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    typedef TImage ImageType;
    typedef typename ImageType::PixelType PixelType;
    typedef typename ImageType::RegionType ImageRegionType;
    typedef PixelType ScalarType;

    itkNewMacro(Self);
    itkTypeMacro(ProlateDistanceToPointImageFilter, ImageToImageFilter);

    itkStaticConstMacro( ImageDimension, unsigned int, ImageType::ImageDimension);

    typedef Vector<ScalarType, 3>              DisplacementType;
    typedef Image<DisplacementType, 3>    DisplacementFieldType;
    typedef ProlateSpheroidalTransform<ScalarType> TransformType;
    typedef typename ImageType::PointType PointType;
    typedef Matrix<ScalarType, 3, 3> BandwidthMatrixType;
    typedef CovariantVector<ScalarType, 3> CovariantVectorType;
    typedef LinearInterpolateImageFunction<DisplacementFieldType, ScalarType> InterpolatorType;
    typedef ContinuousIndex< ScalarType, ImageDimension > ContinuousIndexType;
    
    
    itkSetMacro (Point, PointType);
    itkGetMacro (Point, PointType);
    
    
    void SetInverseDisplacementField (typename DisplacementFieldType::Pointer field)
    {
      m_InverseDisplacementField  = field;
      this->Modified();
    }
    
    void SetTransform (typename TransformType::Pointer forward)
    {
      m_Transform = forward;
      this->Modified();
    }
    
    void SetAlpha (double alpha[3])
    {
      
      // the 3 parameters correspond to the eigen values
      // of the diagonal bandwidth matrix H.
      m_BandwidthMatrix.SetIdentity();
      
      for (unsigned int i=0; i<3; i++)
    	m_BandwidthMatrix[i][i] = alpha[i];
      // the inverse bandwidth matrix
      // is estimated.
      m_SqInverseBandwidthMatrix = m_BandwidthMatrix.GetInverse ();
      m_SqInverseBandwidthMatrix *= m_SqInverseBandwidthMatrix;
      this->Modified();
    }
 
  protected:
    ProlateDistanceToPointImageFilter();
    ~ProlateDistanceToPointImageFilter(){};

    void BeforeThreadedGenerateData(void);
    void AfterThreadedGenerateData(void);
    void ThreadedGenerateData(const ImageRegionType &RegionForThread, int threadId);
    
    void PrintSelf(std::ostream& os, Indent indent) const
    {
      Superclass::PrintSelf(os,indent);
    }

    typename TransformType::Pointer m_Transform;
    typename DisplacementFieldType::Pointer m_InverseDisplacementField;
    typename InterpolatorType::Pointer m_InverseInterpolator;
    
    
    BandwidthMatrixType        m_BandwidthMatrix;
    BandwidthMatrixType        m_SqInverseBandwidthMatrix;
    PointType m_Point;
    
    
  private:

  };


} // end of namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkProlateDistanceToPointImageFilter.txx"
#endif

#endif
