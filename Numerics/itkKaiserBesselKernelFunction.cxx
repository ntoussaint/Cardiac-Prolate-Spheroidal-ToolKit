/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkKaiserBesselKernelFunction.cxx,v $
  Language:  C++
  Date:      $Date: 2008-10-17 01:08:45 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkKaiserBesselKernelFunction.h"

namespace itk
{

KaiserBesselKernelFunction::KaiserBesselKernelFunction()
{
  m_WindowSize = 1.0;
  m_Beta = 4.0;
  m_Resolution = 1000.0;
  m_ScalingFactor = 0.0;
  m_ScaleToUnitIntegral = true;
  
  this->ComputeScalingFactor();
}

KaiserBesselKernelFunction::~KaiserBesselKernelFunction()
{
}

/// Set the window size and re-compute the scaling coefficient.
void KaiserBesselKernelFunction::SetWindowSize(double size)
{
  m_WindowSize = size;
  this->ComputeScalingFactor();
}
  
/// Set the beta parameter and re-compute the scaling coefficient.
void KaiserBesselKernelFunction::SetBeta(double beta)
{
  m_Beta = beta;
  this->ComputeScalingFactor();
}

/// compute the scaling coefficient.
void KaiserBesselKernelFunction::ComputeScalingFactor(void)
{
  double integral = 0.0;
  double step = m_WindowSize / (m_Resolution - 1.0);
  double IO_Beta = this->I0 (m_Beta);
  
  for (double u=0.0; u<=m_WindowSize; u+=step)
  {
    double argument =
      m_Beta *
      std::sqrt ( 1.0  - (2.0 * u / m_WindowSize - 1.0) * (2.0 * u / m_WindowSize - 1.0) );
    integral += step * this->I0 (argument) / IO_Beta;
  }

  m_ScalingFactor = 1.0 / integral;
  
}


/// Returns the modified Bessel function I0(x) for any real x.
double KaiserBesselKernelFunction::I0(double x) const
{
  // Accumulate polynomials in double precision.
  double ax,ans;
  double y; 
  if ((ax=std::fabs(x)) < 3.75)
  {
    //Polynomial fit.
    y=x/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
					 +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  }
  else
  {
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
					  +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
									     +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
														  +y*0.392377e-2))))))));
  }
  return ans;
}  
} // namespace itk
