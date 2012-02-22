/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkKaiserBesselKernelFunction.h,v $
  Language:  C++
  Date:      $Date: 2008-10-17 01:08:45 $
  Version:   $Revision: 1.11 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkKaiserBesselKernelFunction_h
#define __itkKaiserBesselKernelFunction_h

#include "itkKernelFunction.h"
#include <cmath>

namespace itk
{

/** \class KaiserBesselKernelFunction
 * \brief KaiserBessel kernel used for density estimation and nonparameteric
 *  regression.
 *
 * This class encapsulates a KaiserBessel smoothing kernel for
 * density estimation or nonparameteric regression.
 * See documentation for KernelFunction for more details.
 *
 * The Kaiser-Bessel function is a two-parameter family of functions
 * used for digital signal processing, and is defined by the formula
 *
      \f[
        w(k) = 
	\left \{ 
	\begin{matrix}
		\frac{I_0(\beta \sqrt{1 - (2k/N-1)^2})} {I_0(\beta)}
	     & \mbox{if } 0 \leq k \leq N \\  \\
		0 & \mbox{otherwise} \\
	\end{matrix} 
	\right.
      \f]

 * \f$I_0\f$ is the Zero-th order Modified Bessel Function of the first kind.
 *
 * \f$\beta\f$ is a positive real number ( usually around 4.0 )
 * 
 * \f$N\f$ is the size of the support (window).
 * 
 * \sa KernelFunction
 *
 * \ingroup Functions
 */
class ITKCommon_EXPORT KaiserBesselKernelFunction : public KernelFunction
{
public:
  /** Standard class typedefs. */
  typedef KaiserBesselKernelFunction      Self;
  typedef KernelFunction              Superclass;
  typedef SmartPointer<Self>          Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(KaiserBesselKernelFunction, KernelFunction); 

  /** Evaluate the function.

      \f[
        w(k) = 
	\left \{ 
	\begin{matrix}
		\frac{I_0(\beta \sqrt{1 - (2k/N-1)^2})} {I_0(\beta)}
	     & \mbox{if } 0 \leq k \leq N \\  \\
		0 & \mbox{otherwise} \\
	\end{matrix} 
	\right.
      \f]
   */
  inline double Evaluate (const double& u) const
  {
    double my_u = u + m_WindowSize/2.0;
    
    if ( (my_u < 0) || (my_u >m_WindowSize))
      return 0.0;
    
    double argument =
      m_Beta *
      std::sqrt ( 1.0  - (2.0 * my_u / m_WindowSize - 1.0) * (2.0 * my_u / m_WindowSize - 1.0) );

    double value = this->I0 (argument) / this->I0 (m_Beta);
    if (m_ScaleToUnitIntegral) value *= m_ScalingFactor;
    return value;
  }
  
  
  itkGetMacro (WindowSize, double);
  itkGetMacro (Beta, double);
  void SetWindowSize (double size);
  void SetBeta (double beta);

  itkGetMacro (ScaleToUnitIntegral, unsigned int);
  itkSetClampMacro (ScaleToUnitIntegral, unsigned int, 0, 1);
  itkBooleanMacro (ScaleToUnitIntegral);
  
protected:
  KaiserBesselKernelFunction();
  ~KaiserBesselKernelFunction();
  void PrintSelf(std::ostream& os, Indent indent) const
    { Superclass::PrintSelf( os, indent ); }  

  void ComputeScalingFactor(void);
  
  /// Returns the modified Bessel function I0(x) for any real x.
  double I0(double x) const;
  
  double m_WindowSize;
  double m_Beta;
  double m_ScalingFactor;
  double m_Resolution;
  unsigned int m_ScaleToUnitIntegral;
  
private:
  KaiserBesselKernelFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#endif
