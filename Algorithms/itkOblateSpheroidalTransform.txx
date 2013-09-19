/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkOblateSpheroidalTransform.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_OblateSpheroidalTransform_txx_
#define _itk_OblateSpheroidalTransform_txx_

#include "itkOblateSpheroidalTransform.h"
#include <itksys/SystemTools.hxx>

#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <limits>
#include <cassert>

namespace itk
{
  /*
   * Constructor. 
   */
  template <class TPixelType>
  OblateSpheroidalTransform<TPixelType>::OblateSpheroidalTransform()
  {
    this->m_Parameters.SetSize (ParametersDimension);
    this->m_Parameters.Fill (0.0);

    this->m_FixedParameters.SetSize (0);
    
    this->m_Forward = 1;
    
    this->m_Center.Fill (0);
    this->m_LongAxisPoint.Fill (0);
    this->m_ShortAxisPoint.Fill (0);
    
    this->m_Focus1.Fill (0);
    this->m_Focus2.Fill (0);
    this->m_InternalTransform.SetIdentity();
    this->m_InternalTransformInverse.SetIdentity ();

    this->m_LongAxisPoint[1] = std::cosh (1.0);
    this->m_ShortAxisPoint[0] = std::sinh (1.0);

    this->ComputeTransformation();
  }

  /*
   * Destructor. 
   */
  template <class TPixelType>
  OblateSpheroidalTransform<TPixelType>::~OblateSpheroidalTransform()
  {
    
  }
  
  template <class TPixelType>
  void OblateSpheroidalTransform<TPixelType>::ComputeTransformation (void)
  {
    unsigned int count = 0;
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = this->m_Center[i];
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = this->m_LongAxisPoint[i];
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = this->m_ShortAxisPoint[i];
    
    this->m_Axis1 = this->m_ShortAxisPoint  - this->m_Center;
    this->m_Axis2 = this->m_LongAxisPoint - this->m_Center;

    this->m_Lambda1 = this->m_Axis1.GetNorm();
    this->m_Lambda2 = this->m_Axis2.GetNorm();
    this->m_Lambda3 = this->m_Lambda2;
    
    this->m_Axis1.Normalize ();
    this->m_Axis2.Normalize ();
    this->m_Axis3 = CrossProduct (this->m_Axis1, this->m_Axis2);
    // We absolutely need the axes of the spheroid to be orthogonal !
    // By trusting Axis1, we can re-compute axis2 by cross-product
    this->m_Axis2 = CrossProduct (this->m_Axis3, this->m_Axis1);
    
    // A is the rigid transformation matrix from R1 to R0
    for (unsigned int i=0; i<3; i++)
    {
      this->m_InternalTransform[i][0] = this->m_Axis1[i];
      this->m_InternalTransform[i][1] = this->m_Axis2[i];
      this->m_InternalTransform[i][2] = this->m_Axis3[i];
      this->m_InternalTransform[i][3] = this->m_Center[i];
    }
    
    // A-1 is then the transformation matrix from R0 to R1
    this->m_InternalTransformInverse = this->m_InternalTransform.GetInverse();
    
    // The coefficients a b and c correspond to the coefficients of the
    // ellipsoid equation : x2/a2 + y2/b2 + z2/c2 = 1, defined in R1
    this->m_Coefficients[0] =  this->m_Lambda1;
    this->m_Coefficients[1] =  this->m_Lambda2;
    this->m_Coefficients[2] =  this->m_Lambda3;
    // so we have a < b = c ==> oblate spheroid, centered in R1 origin, aligned at x axis, in R1;
    double Eccentricity = std::sqrt (this->m_Lambda2 * this->m_Lambda2 - this->m_Lambda1 * this->m_Lambda1) / this->m_Lambda2;
    this->m_Eccentricity = Eccentricity;
    
    this->m_Focus1 = this->m_Center + Eccentricity * this->m_Lambda2 * this->m_Axis2;
    this->m_Focus2 = this->m_Center - Eccentricity * this->m_Lambda2 * this->m_Axis2;
    
    // The semi foci distance (d in Costa paper 1996)
    this->m_SemiFociDistance = this->m_Focus1.EuclideanDistanceTo (this->m_Focus2) / 2.0;
  }  

  template <class TPixelType>
  typename OblateSpheroidalTransform<TPixelType>::PointType OblateSpheroidalTransform<TPixelType>::ToCartesian(PointType xsi) const
  {

    // xsi has constrains and bounds, lets check them;
    if ( ( (xsi[0] < 0.0) ) ||
	 ( (xsi[1] < - vnl_math::pi / 2.0 ) || (xsi[1] > vnl_math::pi / 2.0) ) ||
	 ( (xsi[2] < 0.0) || (xsi[2] > 2 * vnl_math::pi) ) )
    {
      itkDebugMacro (<<"inconsistent Oblate Coordinate : "<<xsi);
    }
    
    // put the point back into the x-aligned aligned spheroid
    UniformVectorType xs;

    xs[0] = this->m_SemiFociDistance * std::sinh(xsi[0]) * std::sin ( xsi [1]);
    xs[1] = this->m_SemiFociDistance * std::cosh(xsi[0]) * std::cos ( xsi [1]) * std::cos (xsi[2]);
    xs[2] = this->m_SemiFociDistance * std::cosh(xsi[0]) * std::cos ( xsi [1]) * std::sin (xsi[2]);
    xs[3] = 1;
    
    xs = this->m_InternalTransform* xs;

    PointType x;
    for (unsigned int i=0; i<3; i++)
      x[i] = xs[i];
    return x;
  }
  
  template <class TPixelType>
  typename OblateSpheroidalTransform<TPixelType>::PointType OblateSpheroidalTransform<TPixelType>::ToProlate(PointType x) const
  {
    PointType xsi;

    // numerical precision for null tests
    const double epsilon = vcl_numeric_limits<double>::epsilon();

    // put the point back into the x-aligned aligned spheroid,
    // i.e. in referential R1
    UniformVectorType xs;
    for (unsigned int i=0; i<3; i++)
    {
      xs[i]=x[i];
    }
    
    xs[3] = 1;
    xs = this->m_InternalTransformInverse* xs;

    // define a few variables for the polynomial in alpha = cos^2(\xsi_2)
    double A = this->m_SemiFociDistance*this->m_SemiFociDistance;
    double B = xs[1]*xs[1] + xs[2]*xs[2];
    double C = xs[0]*xs[0];
    double gamma1 = (A + B + C - std::sqrt((A+B+C)*(A+B+C) - 4*A*B) ) / (2.0*A);
    
    // by definition gamma1 is positive and should be lower than 1
    // as we'll compute acos(sqrt(gamma1)) below.
    // make sure of this. there could be numerical errors.
    if (gamma1 > (1.0 - epsilon))
      gamma1 = 1.0 - epsilon;
    
    VectorType axis;
    
    // axis2 is the y axis, origin of xsi[2] angle.
    VectorType axis2(0.0);
    axis2[1] = 1;
    
    for (unsigned int i=0; i<3; i++)
      axis[i] = xs[i];
    
    // we project our point to the yz plane by nulling the x component.
    axis[0] = 0;
    // then we normalize it --> we are in the xy unit circle
    // if axis.GetNorm() is null then it means we are directly along the
    // axis of revolution, where xsi[2] is undefined.
    // then we impose xsi[2] = 0 by using axis = axis2.
    if (axis.GetNorm() > epsilon)
      axis.Normalize();
    else
      axis = axis2;
    // then xsi[2] is well defined, apart from cos ambiguity.
    xsi[2] = std::acos (axis[1]);
    // cos ambiguity removed by checking the z sign
    if (axis[2] < 0)
      xsi[2] = 2.0 * vnl_math::pi - xsi[2];
    
    // the real xsi[1] is here,
    xsi[1] = std::acos (std::sqrt (gamma1));
    // cos ambiguity removed by checking the x sign
    if (xs[0] < 0)
      xsi[1] = - xsi[1];
    // the real xsi[0] is here,
    double j1;
    if (gamma1 > epsilon)
    {
      j1 = std::sqrt (B / (A * gamma1));
    }
    else
    {
      // if gamma1 is null, it means we are the axis of revolution,
      // then we have :
      // j1 = std::sqrt (C/A + 1);
      // this equation holds everywhere along the axis BUT between focii,
      // where this train terminates. i.e. prolates coord. are not defined...
      if (std::abs (xs[0]) > this->m_SemiFociDistance)
      {
	j1 = std::sqrt (C/A + 1);
      }
      else
      {
	// we reached the singularity axis between focii, impose xsi[0] = 0
	// by using j1 = 0:
	j1 = 0;
	itkWarningMacro (<<"reach singularity axis between foci\n");
      }
    }

    // the function arccosh is taken as sqrt and logs
    xsi[0] = std::log (j1 + std::sqrt (j1 + 1.0) * std::sqrt (j1 - 1.0));
    
    // xsi has constrains and bounds, let's check them;
    if ( ( (xsi[0] < 0.0) ) ||
    	 ( (xsi[1] < - vnl_math::pi / 2.0) || (xsi[1] > vnl_math::pi / 2.0) ) ||
    	 ( (xsi[2] < 0.0) || (xsi[2] > 2 * vnl_math::pi) ) )
    {
      itkWarningMacro (<<"inconsistent Oblate Coordinate : "<<xsi);
    }
    
    // all singularity situations have been handled, even between focii,
    // no need to test for null xsi[0] now.
    return xsi;
  }


  template <class TPixelType>
  vnl_matrix_fixed<TPixelType,3,3> OblateSpheroidalTransform<TPixelType>::GetJacobianWithRespectToCoordinates(const InputPointType  &x) const
  {
    /// \todo REDO THE COMPUTATION FOR OBLATE
    itkWarningMacro (<< "OblateSpheroidalTransform::GetJacobianWithRespectToCoordinates: CAUTION : This method has to change for the oblate case");

    vnl_matrix_fixed<TPixelType,3,3> jacobian(0.0);
    jacobian[0][0] = jacobian[1][1] = jacobian[2][2] = 1.0;
    
    // error overwhich we consider the system to be inconsistent (non- orthogonal)
    const double epsilon_orthogonality = 0.001;

    // recover the Prolate coordinates point
    PointType xsi = this->m_Forward ? this->TransformPoint (x) : x;

    // error underwhich we display a warning because we reached the singularity
    const double epsilon_singularity = 0.001;
        
    if ( (xsi[0] <= epsilon_singularity) || (xsi[1] <= epsilon_singularity) )
    {
      itkWarningMacro (<<"singularity point : "<<xsi<<"\n"
		       <<"jacobian matrix undefined --> set to identity");
      return jacobian;
    }
    
    UniformVectorType ns, as, ds;
    VectorType n, a, d;
    ns[0] = sinh(xsi[0]) * std::cos(xsi[1]);
    ns[1] = cosh(xsi[0]) * std::sin(xsi[1]) * std::cos(xsi[2]);
    ns[2] = cosh(xsi[0]) * std::sin(xsi[1]) * std::sin(xsi[2]);
    ns[3] = 0;
    
    as[0] = - cosh(xsi[0]) * std::sin(xsi[1]);
    as[1] =   sinh(xsi[0]) * std::cos(xsi[1]) * std::cos(xsi[2]);
    as[2] =   sinh(xsi[0]) * std::cos(xsi[1]) * std::sin(xsi[2]);
    as[3] = 0;
    
    ds[0] = 0;
    ds[1] = - sinh(xsi[0]) * std::sin(xsi[1]) * std::sin(xsi[2]);
    ds[2] = sinh(xsi[0]) * std::sin(xsi[1]) * std::cos(xsi[2]);
    ds[3] = 0;
    
    // the contra-variant basis vectors are in R1 referential
    // they have to be transformed to R0 with this->m_InternalTransform
    // Also, they have to be scaled by the
    // constant \f$ f_{1/2} \f$ (this->m_SemiFociDistance).
    ns = this->m_InternalTransform * (this->m_SemiFociDistance * ns);
    as = this->m_InternalTransform * (this->m_SemiFociDistance * as);
    ds = this->m_InternalTransform * (this->m_SemiFociDistance * ds);
    
    // As we are looking for only re-orientation component of
    // the Jacobian Matrix, we simply multiply the columns by
    // the scale factors:
    double h[3] = {1.0, 1.0, 1.0};
    this->EvaluateScaleFactors (xsi.GetDataPointer(), h);
    
    for (unsigned int i=0; i<3; i++)
    {
      n[i] = ns[i] * h[0];
      a[i] = as[i] * h[1];
      d[i] = ds[i] * h[2];
    }

    // The resulting basis vectors are normals and direct now,
    // with a numerical error < 0.0005 in practice
    
    // fill the output jacobian matrix with column vectors
    for (unsigned int i=0; i<3; i++)
    {
      jacobian[i][0] = n[i];
      jacobian[i][1] = a[i];
      jacobian[i][2] = d[i];
    }

    // a warning is shown is the basis is inconsistent (non-orthogonality) 
    // [jacobian] is a local orthonormal prolate spheroid basis at x;
    // and should be orthogonal and normal basis, let's check that:
    if ( std::abs ( vnl_determinant (jacobian) - 1.0 ) > epsilon_orthogonality ) 
    {
      itkWarningMacro (<<"CAUTION !!! oblate local basis vectors are not orthogonal..."<<std::endl
      		       <<"determinant : "<<vnl_determinant (jacobian)<<std::endl);
    }
    
    // transpose the matrix if we are looking for the inverse transformation
    if (!this->m_Forward)
      jacobian.inplace_transpose();

    // NB : the jacobian matrix is orthonormal but is NOT a rotation matrix.
    // Its determinant is always negative egual to -1 (roto-inversion).
    // This is due to the fact that by construction the basis (g1,g2,g3)
    // is orthogonal but NOT direct. however it does not influence any
    // further process in prolate sph.
    // NB2 : THIS HAS CHANGED. the basis is now DIRECT
    // There was no reason why the basis was actually indirect.
    // A sign error was introduced to $ g_3 $ making the basis indirect,
    // this has been corrected, Note that it did not affect any further results.
    return jacobian;
    
  }
  

  template <class TPixelType>
  void OblateSpheroidalTransform<TPixelType>::EvaluateScaleFactors (double xsi[3], double h[3]) const
  {
    /// \todo REDO THE COMPUTATION FOR OBLATE
    itkWarningMacro (<< "OblateSpheroidalTransform::EvaluateScaleFactors: CAUTION : This method has to change for the oblate case");

    // error underwhich we display a warning because we reached the singularity
    const double epsilon_singularity = 0.001;
    if ( (xsi[0] <= epsilon_singularity) || (xsi[1] <= epsilon_singularity) )
    {
      itkWarningMacro (<<"singularity point : "<<xsi<<"\n"
		       <<"scale factors undefined --> set to null");
      h[0] = h[1] = h[2] = 0.0;
      return;
    }
    
    double d = this->m_SemiFociDistance;
    
    h[0] = 1.0 / ( d * std::sqrt( std::sinh (xsi[0]) * std::sinh (xsi[0]) + std::sin (xsi[1]) * std::sin (xsi[1]) ) );
    h[1] = 1.0 / ( d * std::sqrt( std::sinh (xsi[0]) * std::sinh (xsi[0]) + std::sin (xsi[1]) * std::sin (xsi[1]) ) );
    h[2] = 1.0 / ( d * std::sinh (xsi[0]) * std::sin (xsi[1]) );
  }


  

} // end namespace itk


#endif
