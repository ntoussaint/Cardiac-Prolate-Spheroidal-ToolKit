/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkSliceToVolumePlaneConstrainedMutualInformationCostFunction.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef __itkSliceToVolumePlaneConstrainedMutualInformationCostFunction_h
#define __itkSliceToVolumePlaneConstrainedMutualInformationCostFunction_h


#include "itkSliceToVolumeMutualInformationCostFunction.h"

namespace itk
{
  
  /**
     \class SliceToVolumePlaneConstrainedMutualInformationCostFunction
     \brief This class is a base for the CostFunctions returning a 
     single value
     
     \ingroup ImageRegistration
  */
  template <class TVolumeImage, class TSliceImage>
    class SliceToVolumePlaneConstrainedMutualInformationCostFunction : 
  public SliceToVolumeMutualInformationCostFunction <TVolumeImage, TSliceImage>
{
  public:
  /** Standard class typedefs. */
  typedef SliceToVolumePlaneConstrainedMutualInformationCostFunction     Self;
  typedef SliceToVolumeMutualInformationCostFunction<TVolumeImage, TSliceImage> Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  typedef typename Superclass::PixelType  PixelType;
  typedef typename Superclass::ResamplerType ResamplerType;
  typedef typename Superclass::SliceImageType SliceImageType;
  typedef typename Superclass::SliceRegionType SliceRegionType;
  typedef typename Superclass::SliceSizeType  SliceSizeType;
  
  typedef typename Superclass::VolumeSizeType  VolumeSizeType;
  typedef typename Superclass::PointType PointType;
  typedef typename Superclass::VectorType VectorType;
  typedef typename Superclass::DirectionType DirectionType;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( SliceToVolumePlaneConstrainedMutualInformationCostFunction, SliceToVolumeMutualInformationCostFunction );
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  /**  MeasureType typedef.
   *  It defines a type used to return the cost function value. */
  typedef typename Superclass::MeasureType MeasureType;
  /** DerivativeType typedef.
   *  It defines a type used to return the cost function derivative.  */
  typedef typename Superclass::DerivativeType DerivativeType;
  /**  ParametersType typedef.
   *  It defines a position in the optimization search space. */
  typedef typename Superclass::ParametersType ParametersType;
  /**  ParametersType typedef.
   *  It defines a position in the optimization search space. */
  typedef typename Superclass::ParametersValueType ParametersValueType;
  /** Return the number of parameters required to compute 
   *  this cost function.
   *  This method MUST be overloaded by derived classes. */
  virtual unsigned int GetNumberOfParameters(void) const
  { return 2; }
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

  DirectionType GetDirectionOfConstrain (void)
  { return m_DirectionOfConstrain; }
  void SetDirectionOfConstrain (DirectionType direction)
  { m_DirectionOfConstrain = direction; }
  
protected:
  SliceToVolumePlaneConstrainedMutualInformationCostFunction();
  virtual ~SliceToVolumePlaneConstrainedMutualInformationCostFunction();
  void PrintSelf(std::ostream& os, Indent indent) const;

  DirectionType m_DirectionOfConstrain;
  
private:
  SliceToVolumePlaneConstrainedMutualInformationCostFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk



#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSliceToVolumePlaneConstrainedMutualInformationCostFunction.txx"
#endif

#endif
