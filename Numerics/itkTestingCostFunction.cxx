/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkTestingCostFunction.cxx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itkTestingCostFunction_cxx
#define _itkTestingCostFunction_cxx

#include <itkTestingCostFunction.h>
#include <itkPoint.h>

namespace itk
{

TestingCostFunction
::TestingCostFunction()
{
}

TestingCostFunction
::~TestingCostFunction()
{
}
  
TestingCostFunction::MeasureType
TestingCostFunction
::GetValue (const ParametersType & parameters) const
{
  double norm = 0;
  for (unsigned int i=0; i<2; i++)
    norm += parameters[i] * parameters[i];
  norm = std::sqrt (norm);
  
  return norm;
}

} // end namespace itk
#endif
