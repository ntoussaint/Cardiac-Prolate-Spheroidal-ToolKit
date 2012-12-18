/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkEllipsoidalTransform.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_EllipsoidalTransform_txx_
#define _itk_EllipsoidalTransform_txx_

#include "itkEllipsoidalTransform.h"
#include <itksys/SystemTools.hxx>

#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <limits>
#include <cassert>

#include <vnl/algo/vnl_rnpoly_solve.h>

namespace itk
{
  /*
   * Constructor. 
   */
  template <class TPixelType>
  EllipsoidalTransform<TPixelType>::EllipsoidalTransform() :
    Superclass (OutputSpaceDimension, ParametersDimension)
  {
    this->m_Parameters.SetSize (ParametersDimension);
    this->m_Parameters.Fill (0.0);
    
    this->m_FixedParameters.SetSize (0);
    
    this->m_Jacobian.SetSize (OutputSpaceDimension,InputSpaceDimension);
    this->m_Jacobian.Fill (0.0);
    
    m_Forward = 1;
    
    m_Center.Fill (0);
    m_LongAxisPoint.Fill (0);
    m_ShortAxisPoint1.Fill (0);
    m_ShortAxisPoint2.Fill (0);
    
    m_Focus1.Fill (0);
    m_Focus2.Fill (0);
    m_InternalTransform.SetIdentity();
    m_InternalTransformInverse.SetIdentity ();

    m_LongAxisPoint[0] = std::cosh (1.0);
    m_ShortAxisPoint1[1] = std::sinh (1.0);
    m_ShortAxisPoint2[2] = std::sinh (1.0);

    this->ComputeTransformation();
  }


  /*
   * Destructor. 
   */
  template <class TPixelType>
  EllipsoidalTransform<TPixelType>::~EllipsoidalTransform()
  {
    
  }

  template <class TPixelType>
  void EllipsoidalTransform<TPixelType>::ComputeTransformation (void)
  {
    unsigned int count = 0;
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = m_Center[i];
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = m_LongAxisPoint[i];
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = m_ShortAxisPoint1[i];
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = m_ShortAxisPoint2[i];

    m_Axis1 = m_LongAxisPoint  - m_Center;
    m_Axis2 = m_ShortAxisPoint1 - m_Center;
    m_Axis3 = m_ShortAxisPoint2 - m_Center;

    m_Lambda1 = m_Axis1.GetNorm();
    m_Lambda2 = m_Axis2.GetNorm();
    m_Lambda3 = m_Axis3.GetNorm();

    m_Axis1.Normalize ();
    m_Axis2.Normalize ();
    m_Axis3 = CrossProduct (m_Axis1, m_Axis2);
    // We absolutely need the axes of the spheroid to be orthogonal !
    // By trusting Axis1, we can re-compute axis2 by cross-product
    m_Axis2 = CrossProduct (m_Axis3, m_Axis1);

    // A is the rigid transformation matrix from R1 to R0
    for (unsigned int i=0; i<3; i++)
    {
      m_InternalTransform[i][0] = m_Axis1[i];
      m_InternalTransform[i][1] = m_Axis2[i];
      m_InternalTransform[i][2] = m_Axis3[i];
      m_InternalTransform[i][3] = m_Center[i];
    }
    
    // A-1 is then the transformation matrix from R0 to R1
    m_InternalTransformInverse = m_InternalTransform.GetInverse();

    // The coefficients a b and c correspond to the coefficients of the
    // ellipsoid equation : x2/a2 + y2/b2 + z2/c2 = 1, defined in R1
    m_Coefficients[0] =  m_Lambda1;
    m_Coefficients[1] =  m_Lambda2;
    m_Coefficients[2] =  m_Lambda3;
    // so we have a > b > c ==> ellipsoid, centered in R1 origin, aligned at x axis, in R1;
    double Eccentricity = std::sqrt (m_Lambda1 * m_Lambda1 - m_Lambda2 * m_Lambda2) / m_Lambda1;
    m_Eccentricity = Eccentricity;
    
    m_Focus1 = m_Center + Eccentricity * m_Lambda1 * m_Axis1;
    m_Focus2 = m_Center - Eccentricity * m_Lambda1 * m_Axis1;
    
    // The semi foci distance (d in Costa paper 1996)
    m_SemiFociDistance = m_Focus1.EuclideanDistanceTo (m_Focus2) / 2.0;
    
  }
  
  template <class TPixelType>
  typename EllipsoidalTransform<TPixelType>::OutputPointType EllipsoidalTransform<TPixelType>::TransformPoint(const InputPointType  &x ) const
  {
    OutputPointType ret;
    if (m_Forward)
      ret = this->ToEllipsoidal (x);
    else
      ret = this->ToCartesian (x);
    return ret;
  }

  template <class TPixelType>
  typename EllipsoidalTransform<TPixelType>::VectorType EllipsoidalTransform<TPixelType>::TransformVector(const VectorType &v, const InputPointType &p) const
  {
    VectorType ret;
    if (m_Forward)
      ret = this->ToEllipsoidal (v, p);
    else
      ret = this->ToCartesian (v, p);
    return ret;
  }

  template <class TPixelType>
  void EllipsoidalTransform<TPixelType>::SetParameters (const ParametersType &p)
  {
    this->m_Parameters = p;
    unsigned int count = 0;
    for (unsigned int i=0; i<3; i++) m_Center[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_LongAxisPoint[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_ShortAxisPoint1[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_ShortAxisPoint2[i] = this->m_Parameters[count++];

    this->ComputeTransformation();
  }

  template <class TPixelType>
  void EllipsoidalTransform<TPixelType>::SetParametersByValue (const ParametersType &p)
  {
    this->m_Parameters.SetSize (p.GetSize());
    for (unsigned int i=0; i<p.GetSize(); i++)
      this->m_Parameters[i] = p[i];
    unsigned int count = 0;
    for (unsigned int i=0; i<3; i++) m_Center[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_LongAxisPoint[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_ShortAxisPoint1[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_ShortAxisPoint2[i] = this->m_Parameters[count++];

    this->ComputeTransformation();
  }

  template <class TPixelType>
  typename EllipsoidalTransform<TPixelType>::JacobianType & EllipsoidalTransform<TPixelType>::GetJacobian(const InputPointType  &x) const
  {
    // jacobian according to the "parameters" is definitely NOT defined...
    itkWarningMacro (<<"jacobian according to the \"parameters\" is definitely NOT defined..."<<"\n");
    
    return this->m_Jacobian;
  }

  template <class TPixelType>
  vnl_matrix_fixed<TPixelType,3,3> EllipsoidalTransform<TPixelType>::GetJacobianWithRespectToCoordinates(const InputPointType  &x) const
  {
    vnl_matrix_fixed<TPixelType,3,3> jacobian(0.0);
    jacobian[0][0] = jacobian[1][1] = jacobian[2][2] = 1.0;
    
    // error overwhich we consider the system to be inconsistent (non- orthogonal)
    const double epsilon_orthogonality = 0.001;
    
    // recover the Ellipsoidal coordinates point
    OutputPointType xsi = m_Forward ? this->TransformPoint (x) : x;
    
    const double epsilon_singularity = 0.0001;
    // error underwhich we display a warning because we reached the singularity
    if ( ( std::abs (xsi[0] - xsi[1]) < epsilon_singularity ) ||
	 ( std::abs (xsi[1] - xsi[2]) < epsilon_singularity ) ||
	 ( std::abs (xsi[2] - xsi[0]) < epsilon_singularity ) )
    {
      itkWarningMacro (<<"singularity point : "<<xsi<<"\n"
		       <<"jacobian undefined --> set to Id");
    }
    
    double a2 = m_Lambda1*m_Lambda1;
    double b2 = m_Lambda2*m_Lambda2;
    double c2 = m_Lambda3*m_Lambda3;

    double A = (a2 - b2)*(a2 - c2);
    double B = (b2 - c2)*(b2 - a2);
    double C = (c2 - a2)*(c2 - b2);

    double B11 = ( (a2 - xsi[1]) * (a2 - xsi[2]) ) / A;
    double B12 = ( (a2 - xsi[2]) * (a2 - xsi[0]) ) / A;
    double B13 = ( (a2 - xsi[0]) * (a2 - xsi[1]) ) / A;

    double B21 = ( (b2 - xsi[1]) * (b2 - xsi[2]) ) / B;
    double B22 = ( (b2 - xsi[2]) * (b2 - xsi[0]) ) / B;
    double B23 = ( (b2 - xsi[0]) * (b2 - xsi[1]) ) / B;
    
    double B31 = ( (c2 - xsi[1]) * (c2 - xsi[2]) ) / C;
    double B32 = ( (c2 - xsi[2]) * (c2 - xsi[0]) ) / C;
    double B33 = ( (c2 - xsi[0]) * (c2 - xsi[1]) ) / C;
    
    UniformVectorType ns, as, ds;
    VectorType n, a, d;

    // need to choose the signs with xsi[3]...
    
    ns[0] = std::sqrt( B11 / (4.0 * (a2 - xsi[0])) );
    ns[0] = std::sqrt( B12 / (4.0 * (b2 - xsi[1])) );
    ns[0] = std::sqrt( B13 / (4.0 * (c2 - xsi[2])) );
    ns[3] = 0;
    
    as[0] = std::sqrt( B21 / (4.0 * (a2 - xsi[0])) );
    as[0] = std::sqrt( B22 / (4.0 * (b2 - xsi[1])) );
    as[0] = std::sqrt( B23 / (4.0 * (c2 - xsi[2])) );
    as[3] = 0;
    
    ds[0] = std::sqrt( B31 / (4.0 * (a2 - xsi[0])) );
    ds[0] = std::sqrt( B32 / (4.0 * (b2 - xsi[1])) );
    ds[0] = std::sqrt( B33 / (4.0 * (c2 - xsi[2])) );
    ds[3] = 0;
    
    // the contra-variant basis vectors are in R1 referential
    // they have to be transformed to R0 with m_InternalTransform
    ns = m_InternalTransform * ns;
    as = m_InternalTransform * as;
    ds = m_InternalTransform * ds;
    
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
      itkWarningMacro (<<"CAUTION !!! prolate local basis vectors are not orthogonal..."<<std::endl
      		       <<"determinant : "<<vnl_determinant (jacobian)<<std::endl);
    }
    
    // transpose the matrix if we are looking for the inverse transformation
    if (!m_Forward)
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
  void EllipsoidalTransform<TPixelType>::EvaluateLocalBasis(InputPointType x,
							    VectorType &n,
							    VectorType &a,
							    VectorType &d) const
  {
    vnl_matrix_fixed<double,3,3> jacobian = this->GetJacobianWithRespectToCoordinates (x);
    // The basis vectors correspond to each column of the jacobian matrix
    for (unsigned int i=0; i<3; i++)
    {
      n[i] = jacobian[i][0];
      a[i] = jacobian[i][1];
      d[i] = jacobian[i][2];
    }
  }
  
  template <class TPixelType>
  bool EllipsoidalTransform<TPixelType>::GetInverse( Self* inverse) const
  {
    inverse->SetParameters (this->m_Parameters);
    inverse->SetForward (!this->GetForward());
    return true;
  }
  
  template <class TPixelType>
  void EllipsoidalTransform<TPixelType>::SetIdentity(void)
  {
    m_Center.Fill (0.0);
    m_LongAxisPoint.Fill (0.0);
    m_LongAxisPoint[0] = std::cosh (1.0);
    m_ShortAxisPoint1.Fill (0.0);
    m_ShortAxisPoint2.Fill (0.0);
    m_ShortAxisPoint1[1] = std::sinh (1.0);
    m_ShortAxisPoint2[2] = std::sinh (1.0);

    this->ComputeTransformation();
  }

  template <class TPixelType>
  typename EllipsoidalTransform<TPixelType>::OutputPointType EllipsoidalTransform<TPixelType>::ToCartesian(InputPointType xsi) const
  {
    double a2 = m_Lambda1*m_Lambda1;
    double b2 = m_Lambda2*m_Lambda2;
    double c2 = m_Lambda3*m_Lambda3;
    // xsi has bounds, lets check them;
    if ( ( (xsi[0] > c2) ) ||
	 ( (xsi[1] < c2) || (xsi[1] > b2) ) ||
	 ( (xsi[2] < b2) || (xsi[2] > a2) ) )
    {
      itkWarningMacro (<<"inconsistent Ellipsoidal Coordinate : "<<xsi);
    }
    
    // put the point back into the x-aligned aligned spheroid
    UniformVectorType xs;
    ScalarType a = m_Lambda1, b = m_Lambda2, c = m_Lambda3;    
    
    int sign = (int)xsi[3];
    
    double e2 = (sign >= 4) ? 1 : -1;
    if (e2>0) sign -= 4;
    double e1 = (sign >= 2) ? 1 : -1;
    double e0 = ( ((int)sign % 2) != 0 ) ? 1 : -1;
    
    // need to choose the sign with xsi[3]
    xs[0] = e0 * std::sqrt ( ( (a*a - xsi[0]) * (a*a - xsi[1]) * (a*a - xsi[2]) ) / ( (a*a - b*b) * (a*a - c*c) ) );
    xs[1] = e1 * std::sqrt ( ( (b*b - xsi[0]) * (b*b - xsi[1]) * (b*b - xsi[2]) ) / ( (b*b - a*a) * (b*b - c*c) ) );
    xs[2] = e2 * std::sqrt ( ( (c*c - xsi[0]) * (c*c - xsi[1]) * (c*c - xsi[2]) ) / ( (c*c - a*a) * (c*c - b*b) ) );
    
    xs[3] = 1;
    
    xs = m_InternalTransform* xs;
    
    OutputPointType x;
    for (unsigned int i=0; i<3; i++)
      x[i] = xs[i];

    x[3] = xsi[3];
    
    return x;
  }
  
  template <class TPixelType>
  typename EllipsoidalTransform<TPixelType>::OutputPointType EllipsoidalTransform<TPixelType>::ToEllipsoidal(InputPointType x) const
  {
    const double epsilon = 0.01;
    
    OutputPointType xsi; xsi[0] = xsi[1] = xsi[2] = 0.0; xsi[3] = 0.0;
    
    // put the point back into the x-aligned aligned spheroid,
    // i.e. in referential R1
    UniformVectorType xs;
    for (unsigned int i=0; i<3; i++)
      xs[i]=x[i];
    xs[3] = 1;
    xs = m_InternalTransformInverse* xs;
    
    // define variables for the polynomial in \xsi
    double a=m_Lambda1, b=m_Lambda2, c=m_Lambda3;
    if (std::abs (b - a) < epsilon)
      b = a - epsilon;
    if (std::abs (c - b) < epsilon)
      c = b - epsilon;
    
    double a2 = a*a, b2=b*b, c2=c*c;
    double x2 = xs[0]*xs[0], y2 = xs[1]*xs[1], z2 = xs[2]*xs[2];
    double A = a2 + b2 + c2;
    double B = a2 * b2 + b2 * c2 + a2 * c2;
    double C = a2 * b2 * c2;
    double D = x2 + y2 + z2;
    double E = x2 * (b2 + c2) + y2 * (a2 + c2) + z2 * (a2 + b2);
    double F = x2 *  b2 * c2  + y2 *  a2 * c2  + z2 *  a2 * b2 ;
    
    double s2 = D-A;
    double s1 = B-E;
    double s0 = F-C;
    double Q = (3.0*s1 - s2*s2) / 9.0;
    double R = (9.0*s2*s1 - 27.0*s0 - 2.0*s2*s2*s2) / 54.0;
    double delta = Q*Q*Q + R*R;
    
    if (delta > 0)
    {
      itkWarningMacro (<<"found non-real roots at  xs = "<<xs<<" with Delta = "<<delta);
      return xsi;
    }

    double theta = std::acos( R / std::sqrt (-Q*Q*Q) );
    
    std::vector<ScalarType> sorter;
    sorter.push_back (2.0 * std::sqrt (-Q) * std::cos( (theta + 0.0 * vnl_math::pi) / 3.0 ) - s2 / 3.0);
    sorter.push_back (2.0 * std::sqrt (-Q) * std::cos( (theta + 2.0 * vnl_math::pi) / 3.0 ) - s2 / 3.0);
    sorter.push_back (2.0 * std::sqrt (-Q) * std::cos( (theta + 4.0 * vnl_math::pi) / 3.0 ) - s2 / 3.0);
    
    std::sort (sorter.begin(), sorter.end());
    for (unsigned int i=0; i<sorter.size(); i++)
      xsi[i] = sorter[i];

    // put xsi[3] according to the signs of xs[0], xs[1], xs[2]
    xsi[3] =
      (double)(xs[0] >= 0) * std::pow (2,0) +
      (double)(xs[1] >= 0) * std::pow (2,1) +
      (double)(xs[2] >= 0) * std::pow (2,2);
    
    return xsi;
  }

  
  template <class TPixelType>
  typename EllipsoidalTransform<TPixelType>::VectorType EllipsoidalTransform<TPixelType>::ToCartesian(VectorType v, InputPointType p) const
  { 
    double a2 = m_Lambda1*m_Lambda1;
    double b2 = m_Lambda2*m_Lambda2;
    double c2 = m_Lambda3*m_Lambda3;
    
    // xsi has bounds, lets check them;
    if ( ( (p[0] > c2) ) ||
	 ( (p[1] < c2) || (p[1] > b2) ) ||
	 ( (p[2] < b2) || (p[2] > a2) ) )
    {
      itkWarningMacro (<<"inconsistent Ellipsoidal Coordinate : "<<p);
    }
    
    MatrixType matrix; matrix = this->GetJacobianWithRespectToCoordinates(p);
    VectorType ret = matrix * v;
    
    return ret;
  }

  template <class TPixelType>
  typename EllipsoidalTransform<TPixelType>::VectorType EllipsoidalTransform<TPixelType>::ToEllipsoidal(VectorType v, InputPointType p) const
  {
    MatrixType matrix; matrix = this->GetJacobianWithRespectToCoordinates(p);
    MatrixType transpose; transpose = matrix.GetTranspose ();
    VectorType ret = transpose * v;
    
    return ret;
  }

  template <class TPixelType>
  void EllipsoidalTransform<TPixelType>::EvaluateScaleFactors (double xsi[3], double h[3]) const
  {
    
    const double epsilon_singularity = 0.0001;
    // error underwhich we display a warning because we reached the singularity
    if ( ( std::abs (xsi[0] - xsi[1]) < epsilon_singularity ) ||
	 ( std::abs (xsi[1] - xsi[2]) < epsilon_singularity ) ||
	 ( std::abs (xsi[2] - xsi[0]) < epsilon_singularity ) )
    {
      itkWarningMacro (<<"singularity point : "<<xsi<<"\n"
		       <<"scale factors undefined --> set to null");
      h[0] = h[1] = h[2] = 0.0;
      return;
    }

    double a2 = m_Lambda1*m_Lambda1;
    double b2 = m_Lambda2*m_Lambda2;
    double c2 = m_Lambda3*m_Lambda3;
    
    h[0] = std::sqrt( ( 4.0 * (a2 - xsi[0]) * (b2 - xsi[0]) * (c2 - xsi[0]) ) / ( (xsi[1] - xsi[0]) * (xsi[2] - xsi[0]) ) );
    h[1] = std::sqrt( ( 4.0 * (a2 - xsi[1]) * (b2 - xsi[1]) * (c2 - xsi[1]) ) / ( (xsi[2] - xsi[1]) * (xsi[0] - xsi[1]) ) );
    h[2] = std::sqrt( ( 4.0 * (a2 - xsi[2]) * (b2 - xsi[2]) * (c2 - xsi[2]) ) / ( (xsi[0] - xsi[2]) * (xsi[1] - xsi[2]) ) );
  }

  template <class TPixelType>
  void EllipsoidalTransform<TPixelType>::PrintSelf(std::ostream& os, Indent indent) const
  {
    this->Superclass::PrintSelf(os,indent);

    os << indent << "Excentricity : "<<m_Eccentricity << "\n";
    os << indent << "Semi Foci Distance : "<<m_SemiFociDistance << "\n";
    
    os << indent << "Coefficients: " 
       << "\n\ta: " << m_Coefficients[0]
       << "\n\tb: " << m_Coefficients[1]
       << "\n\tc: " << m_Coefficients[2] << "\n";
    os << indent << "Lambdas: " 
       << "\n\tl1: " << m_Lambda1
       << "\n\tl2: " << m_Lambda2
       << "\n\tl3: " << m_Lambda3 << "\n";
    os << indent << "Center: " 
       << "\n\tO1: " << m_Center[0]
       << "\n\tO2: " << m_Center[1]
       << "\n\tO3: " << m_Center[2] << "\n";
    os << indent << "v1: " 
       << "\n\tv11: " << m_Axis1[0]
       << "\n\tv12: " << m_Axis1[1]
       << "\n\tv13: " << m_Axis1[2] << "\n";
    os << indent << "v2: " 
       << "\n\tv21: " << m_Axis2[0]
       << "\n\tv22: " << m_Axis2[1]
       << "\n\tv23: " << m_Axis2[2] << "\n";
    os << indent << "v3: " 
       << "\n\tv31: " << m_Axis3[0]
       << "\n\tv32: " << m_Axis3[1]
       << "\n\tv33: " << m_Axis3[2] << "\n";
  }

  /*
   * Evaluate ellipsoid equation.
   */
  template <class TPixelType>
  double EllipsoidalTransform<TPixelType>::EvaluateFunction(PointType x) const
  {
    // Transform input point from R0 to the R1 scheme,
    // where the ellipsoid equation is defined
    UniformVectorType xs;
    for (unsigned int i=0; i<3; i++)
      xs[i]=x[i];
    xs[3] = 1;

    xs = m_InternalTransformInverse* xs;

    const double *a = m_Coefficients;
  
    // ellipsoid formula : x2/a2 + y2/b2 + z2/c2 = 1
    return ( xs[0]*xs[0] / (a[0]*a[0]) + xs[1]*xs[1] / (a[1]*a[1]) + xs[2]*xs[2] / (a[2]*a[2]) - 1 );
  }


  /*
   * Evaluate ellipsoid equation.
   */
  template <class TPixelType>
  double EllipsoidalTransform<TPixelType>::EvaluateFunction(double x[3]) const
  {
    PointType xi;
    for (unsigned int i=0; i<3; i++)
      xi[i] = x[i];
    return this->EvaluateFunction (xi);
  }
  
  /*
   * Evaluate ellipsoid equation gradient.
   */
  template <class TPixelType>
  void EllipsoidalTransform<TPixelType>::EvaluateGradient(PointType x, VectorType &n) const
  {

    // Transform input point from R0 to the R1 frame,
    // where the ellipsoid equation is defined
    UniformVectorType xs;
    for (unsigned int i=0; i<3; i++)
      xs[i]=x[i];
    xs[3] = 1;

    xs = m_InternalTransformInverse * xs;

    const double *a=m_Coefficients;
    // define the gradient at xs (in R1)
    UniformVectorType ns;

    for (unsigned int i=0; i<3; i++)
      ns[i] = 2.0 * ( 1.0 / (a[i]*a[i]) ) * xs[i];

    // ns is a vector
    ns[3] = 0;
    ns = m_InternalTransform * ns;
    
    n[0] = ns[0];
    n[1] = ns[1];
    n[2] = ns[2];
  }


  template <class TPixelType>
  double EllipsoidalTransform<TPixelType>::EstimateGeodesicLength1(InputPointType xi, InputVectorType dxi, unsigned int divisions) const
  {
    InputVectorType step = dxi / (double)(divisions);
    double length = 0.0;
    InputPointType point1 = xi, point2 = xi;
    InputPointType cart1, cart2;
    InputVectorType dL;
    
    for (unsigned int i=0; i<divisions; i++)
    {
      point2 = point1 + step;
      cart1 = this->ToCartesian (point1);
      cart2 = this->ToCartesian (point2);
      dL = cart2 - cart1;
      length += dL.GetNorm();
      point1 = point2;
    }

    return length;
  }

  template <class TPixelType>
  double EllipsoidalTransform<TPixelType>::EstimateGeodesicLength2(InputPointType xi, InputVectorType dxi, unsigned int divisions) const
  {
    InputVectorType step = dxi / (double)(divisions);
    InputPointType pt = xi;
    double h[3];
    
    double length = 0.0;
    for (unsigned int i=0; i<divisions; i++)
    {
      this->EvaluateScaleFactors (pt.GetDataPointer(), h);
      double ds = 0.0;
      for (unsigned int j=0; j<3; j++) ds += ( step[j] * step[j] ) / ( h[j] * h[j] );
      ds = std::sqrt ( ds );
      length += ds;
      pt += step;
    }
    
    return length;
  }


} // end namespace itk


#endif
