/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkSliceToVolumeMutualInformationCostFunction.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef __itkSliceToVolumeMutualInformationCostFunction_h
#define __itkSliceToVolumeMutualInformationCostFunction_h

#include "itkSingleValuedCostFunction.h"
#include "itkArray.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"
#include "itkMatrix.h"
#include "itkImage.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkTranslationTransform.h"
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkCorrelationCoefficientHistogramImageToImageMetric.h>
#include <vector>
namespace itk
{
  
  /**
     \class SliceToVolumeMutualInformationCostFunction
     \brief This class is a base for the CostFunctions returning a 
     single value
     
     \ingroup ImageRegistration
  */
  template <class TVolumeImage, class TSliceImage>
    class ITK_EXPORT SliceToVolumeMutualInformationCostFunction : 
  public SingleValuedCostFunction 
{
public:
  /** Standard class typedefs. */
  typedef SliceToVolumeMutualInformationCostFunction     Self;
  typedef SingleValuedCostFunction            Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;
  typedef TSliceImage                         SliceImageType;
  typedef typename SliceImageType::Pointer    SliceImagePointerType;
  typedef TVolumeImage                        VolumeImageType;
  typedef typename VolumeImageType::Pointer   VolumeImagePointerType;
  typedef typename VolumeImageType::PixelType  PixelType;
  typedef ResampleImageFilter<VolumeImageType, VolumeImageType> ResamplerType;
  typedef typename ResamplerType::Pointer ResamplerPointerType;
  typedef LinearInterpolateImageFunction<VolumeImageType, PixelType> VolumeInterpolatorType;
  typedef typename VolumeInterpolatorType::Pointer VolumeInterpolatorPointerType;
  typedef LinearInterpolateImageFunction<SliceImageType, PixelType> SliceInterpolatorType;
  typedef typename SliceInterpolatorType::Pointer SliceInterpolatorPointerType;
  typedef typename VolumeImageType::PixelContainer PixelContainerType;
  typedef typename SliceImageType::RegionType SliceRegionType;
  typedef typename SliceImageType::SizeType  SliceSizeType;
  typedef typename VolumeImageType::SizeType  VolumeSizeType;
  typedef typename VolumeImageType::PointType PointType;
  
  typedef TranslationTransform<double, 3> TransformType;
  typedef TransformType::Pointer TransformPointerType;
  typedef TranslationTransform<double, 2> SliceTransformType;
  typedef SliceTransformType::Pointer SliceTransformPointerType;

  typedef itk::Vector<PixelType, 3> VectorType;
  typedef typename VolumeImageType::DirectionType DirectionType;
  
  typedef MattesMutualInformationImageToImageMetric<SliceImageType, SliceImageType> MetricType;
  typedef typename MetricType::Pointer MetricPointerType;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( SliceToVolumeMutualInformationCostFunction, SingleValuedCostFunction );
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  /**  MeasureType typedef.
   *  It defines a type used to return the cost function value. */
  typedef Superclass::MeasureType MeasureType;
  /** DerivativeType typedef.
   *  It defines a type used to return the cost function derivative.  */
  typedef Superclass::DerivativeType DerivativeType;
  /**  ParametersType typedef.
   *  It defines a position in the optimization search space. */
  typedef Superclass::ParametersType ParametersType;
  /**  ParametersType typedef.
   *  It defines a position in the optimization search space. */
  typedef Superclass::ParametersValueType ParametersValueType;
  /** Return the number of parameters required to compute 
   *  this cost function.
   *  This method MUST be overloaded by derived classes. */
  virtual unsigned int GetNumberOfParameters(void) const
  { return 3; }
  /** This method returns the derivative of the cost function corresponding
    * to the specified parameters.   */ 
  virtual void GetDerivative( const ParametersType & parameters,
                              DerivativeType & derivative ) const
  {
    ///\todo Derive this criterion...
  }
  /** This method returns the value of the cost function corresponding
    * to the specified parameters.    */ 
  virtual MeasureType GetValue( const ParametersType & parameters ) const;

  /** Connect the Volume Image.  */
  itkSetObjectMacro( VolumeImage, VolumeImageType );

  /** Get the Volume Image. */
  itkGetObjectMacro( VolumeImage, VolumeImageType );

  /** Connect the Slice Image.  */
  itkSetObjectMacro( SliceImage, VolumeImageType );

  /** Get the Slice Image. */
  itkGetObjectMacro( SliceImage, VolumeImageType );
  
protected:
  SliceToVolumeMutualInformationCostFunction();
  virtual ~SliceToVolumeMutualInformationCostFunction();
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  VolumeInterpolatorPointerType m_VolumeInterpolator;
  SliceInterpolatorPointerType m_SliceInterpolator;
  ResamplerPointerType m_Resampler;
  MetricPointerType m_Metric;
  VolumeImagePointerType m_VolumeImage;
  VolumeImagePointerType m_SliceImage;
  SliceImagePointerType m_Input1;
  SliceImagePointerType m_Input2;
  TransformPointerType m_Transform;
  TransformPointerType m_InverseTransform;
  SliceTransformPointerType m_NullTransform;
  unsigned int m_ParametersDimension;
  
private:
  SliceToVolumeMutualInformationCostFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk



#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSliceToVolumeMutualInformationCostFunction.txx"
#endif

#endif
