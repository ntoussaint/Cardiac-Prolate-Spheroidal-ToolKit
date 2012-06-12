/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkTensorMeshImageHybridCostFunction.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef __itkTestingCostFunction_h
#define __itkTestingCostFunction_h

#include <itkSingleValuedCostFunction.h>

namespace itk
{
  
class ITK_EXPORT TestingCostFunction : 
    public SingleValuedCostFunction 
{
public:
  /** Standard class typedefs. */
  typedef TestingCostFunction  Self;
  typedef SingleValuedCostFunction            Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;
  itkTypeMacro( TestingCostFunction, SingleValuedCostFunction );
  itkNewMacro(Self);

  typedef Superclass::MeasureType MeasureType;
  typedef Superclass::DerivativeType DerivativeType;
  typedef Superclass::ParametersType ParametersType;
  typedef Superclass::ParametersValueType ParametersValueType;

  typedef double ScalarType;
  virtual unsigned int GetNumberOfParameters(void) const
  { return 3; }
  
  virtual MeasureType GetValue      ( const ParametersType & parameters ) const;
  virtual void        GetDerivative ( const ParametersType & parameters, DerivativeType & derivative ) const
  { itkExceptionMacro (<<"NOT IMPLEMENTED"); }
  
protected:
  TestingCostFunction();
  virtual ~TestingCostFunction();
  
private:
  TestingCostFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
};

} // end namespace itk



#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTestingCostFunction.cxx"
#endif

#endif
