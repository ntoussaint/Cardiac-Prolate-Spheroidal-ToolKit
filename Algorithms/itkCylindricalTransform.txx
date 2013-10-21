/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkCylindricalTransform.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_CylindricalTransform_txx_
#define _itk_CylindricalTransform_txx_

#include "itkCylindricalTransform.h"
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
  CylindricalTransform<TPixelType>::CylindricalTransform() :
    Superclass (OutputSpaceDimension)
  {
    this->m_Parameters.SetSize (ParametersDimension);
    this->m_Parameters.Fill (0.0);

    this->m_FixedParameters.SetSize (0);
    
    m_Forward = 1;
    
    m_Center.Fill (0);
    m_LongAxisPoint.Fill (0);
    m_ShortAxisPoint.Fill (0);
    
    m_InternalTransform.SetIdentity();
    m_InternalTransformInverse.SetIdentity ();

    m_LongAxisPoint[0] = std::cosh (1.0);
    m_ShortAxisPoint[1] = std::sinh (1.0);

    this->ComputeTransformation();
  }


  /*
   * Destructor. 
   */
  template <class TPixelType>
  CylindricalTransform<TPixelType>::~CylindricalTransform()
  {
    
  }


  template <class TPixelType>
  void CylindricalTransform<TPixelType>::ComputeTransformation (void)
  {
    unsigned int count = 0;
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = m_Center[i];
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = m_LongAxisPoint[i];
    for (unsigned int i=0; i<3; i++) this->m_Parameters[count++] = m_ShortAxisPoint[i];

    m_Axis1 = m_LongAxisPoint  - m_Center;
    m_Axis2 = m_ShortAxisPoint - m_Center;

    m_Axis1.Normalize ();
    m_Axis2.Normalize ();
    m_Axis3 = CrossProduct (m_Axis1, m_Axis2);
    // We absolutely need the axes of the cylindrical coord. to be orthogonal !
    // By trusting Axis1, we can re-compute axis2 by cross-product
    m_Axis2 = CrossProduct (m_Axis3, m_Axis1);

    // A is the rigid transformation matrix from R1 to R0
    for (unsigned int i=0; i<3; i++)
    {
      m_InternalTransform[i][0] = m_Axis2[i];
      m_InternalTransform[i][1] = m_Axis3[i];
      m_InternalTransform[i][2] = m_Axis1[i];
      m_InternalTransform[i][3] = m_Center[i];
    }
    // A-1 is then the transformation matrix from R0 to R1
    m_InternalTransformInverse = m_InternalTransform.GetInverse();

    // R1 is a referential centered in m_Center and with an x-axis towards the short axis point,
    // the y-axis on the short axis plane, and the z-axis along the ventricle axis, towards the apex.    

  }
  
  template <class TPixelType>
  typename CylindricalTransform<TPixelType>::OutputPointType CylindricalTransform<TPixelType>::TransformPoint(const InputPointType  &x ) const
  {
    OutputPointType ret;
    if (m_Forward)
      ret = this->ToCylindrical (x);
    else
      ret = this->ToCartesian (x);
    return ret;
  }

  template <class TPixelType>
  typename CylindricalTransform<TPixelType>::OutputVectorType CylindricalTransform<TPixelType>::TransformVector(const InputVectorType &v, const InputPointType &p) const
  {
    OutputVectorType ret;
    if (m_Forward)
      ret = this->ToCylindrical (v, p);
    else
      ret = this->ToCartesian (v, p);
    return ret;
  }

  template <class TPixelType>
  typename CylindricalTransform<TPixelType>::OutputVnlVectorType CylindricalTransform<TPixelType>::TransformVector(const InputVnlVectorType &v, const InputPointType &p) const
  {
    VectorType mv, mret;
    OutputVnlVectorType ret;
    for (unsigned int i=0; i<3; i++) mv[i] = v[i];
    mret = this->TransformVector (mv, p);
    for (unsigned int i=0; i<3; i++) ret[i] = mret[i];
    return ret;
  }

  template <class TPixelType>
  typename CylindricalTransform<TPixelType>::OutputCovariantVectorType CylindricalTransform<TPixelType>::TransformCovariantVector(const InputCovariantVectorType &v, const InputPointType &p) const
  {
    VectorType mv, mret;
    OutputCovariantVectorType ret;
    for (unsigned int i=0; i<3; i++) mv[i] = v[i];
    mret = this->TransformVector (mv, p);
    for (unsigned int i=0; i<3; i++) ret[i] = mret[i];
    return ret;
  }

  template <class TPixelType>
  void CylindricalTransform<TPixelType>::SetParameters (const ParametersType &p)
  {
    this->m_Parameters = p;
    unsigned int count = 0;
    for (unsigned int i=0; i<3; i++) m_Center[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_LongAxisPoint[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_ShortAxisPoint[i] = this->m_Parameters[count++];

    this->ComputeTransformation();
  }

  template <class TPixelType>
  void CylindricalTransform<TPixelType>::SetParametersByValue (const ParametersType &p)
  {
    this->m_Parameters.SetSize (p.GetSize());
    for (unsigned int i=0; i<p.GetSize(); i++)
      this->m_Parameters[i] = p[i];
    unsigned int count = 0;
    for (unsigned int i=0; i<3; i++) m_Center[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_LongAxisPoint[i] = this->m_Parameters[count++];
    for (unsigned int i=0; i<3; i++) m_ShortAxisPoint[i] = this->m_Parameters[count++];

    this->ComputeTransformation();
  }

  template <class TPixelType>
  vnl_matrix_fixed<TPixelType,3,3> CylindricalTransform<TPixelType>::GetJacobianWithRespectToCoordinates(const InputPointType  &x) const
  {

    vnl_matrix_fixed<TPixelType,3,3> jacobian(0.0);
    jacobian[0][0] = jacobian[1][1] = jacobian[2][2] = 1.0;
    
    // error overwhich we consider the system to be inconsistent (non- orthogonal)
    const double epsilon_orthogonality = 0.001;

    // recover the Cylindrical coordinates point
    PointType xsi = m_Forward ? this->TransformPoint (x) : x;

    // error underwhich we display a warning because we reached the singularity
    const double epsilon_singularity = 0.001;
    
    if ( (xsi[0] <= epsilon_singularity) )
    {
      itkWarningMacro (<<"singularity point : "<<xsi<<"\n"
		       <<"jacobian matrix undefined --> set to identity");
      return jacobian;
    }
    
    UniformVectorType ns, as, ds;
    VectorType n, a, d;
    // dx / dxsi[0]
    ns[0] = cos (xsi[1]);
    // dy / dxsi[0]
    ns[1] = sin (xsi[1]);
    // dz / dxsi[0]
    ns[2] = 0;
    ns[3] = 0;
    
    // dx / dxsi[1]
    as[0] = - xsi[0] * sin (xsi[1]);
    // dy / dxsi[1]
    as[1] =   xsi[0] * cos (xsi[1]);
    // dz / dxsi[1]
    as[2] =   0;
    as[3] = 0;
    
    // dx / dxsi[2]
    ds[0] = 0;
    // dy / dxsi[2]
    ds[1] = 0;
    // dz / dxsi[2]
    ds[2] = 1;
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
  void CylindricalTransform<TPixelType>::EvaluateLocalBasis(InputPointType x,
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
  bool CylindricalTransform<TPixelType>::GetInverse( Self* inverse) const
  {
    inverse->SetParameters (this->m_Parameters);
    inverse->SetForward (!this->GetForward());
    return true;
  }
  
  template <class TPixelType>
  void CylindricalTransform<TPixelType>::SetIdentity(void)
  {
    m_Center.Fill (0.0);
    m_LongAxisPoint.Fill (0.0);
    m_ShortAxisPoint.Fill (0.0);
    m_LongAxisPoint[0] = std::cosh (1.0);
    m_ShortAxisPoint[1] = std::sinh (1.0);

    this->ComputeTransformation();
  }

  template <class TPixelType>
  typename CylindricalTransform<TPixelType>::PointType CylindricalTransform<TPixelType>::ToCartesian(PointType xsi) const
  {

    
    // xsi has constrains and bounds, lets check them;
    if ( ( (xsi[0] < 0.0) ) ||
	 ( (xsi[1] < - vnl_math::pi / 2.0 ) || (xsi[1] > vnl_math::pi / 2.0) ) )
    {
      itkDebugMacro (<<"inconsistent Cylindrical Coordinate : "<<xsi);
    }
    
    // put the point back into the x-aligned aligned spheroid
    UniformVectorType xs;

    xs[0] = xsi[0] * cos (xsi[1]);
    xs[1] = xsi[0] * sin (xsi[1]);
    xs[2] = xsi[2];
    xs[3] = 1;
    
    xs = m_InternalTransform * xs;

    PointType x;
    for (unsigned int i=0; i<3; i++)
      x[i] = xs[i];
    return x;
  }


  
  template <class TPixelType>
  typename CylindricalTransform<TPixelType>::PointType CylindricalTransform<TPixelType>::ToCylindrical(PointType x) const
  {
    PointType xsi;

    // numerical precision for null tests
    const double epsilon = vcl_numeric_limits<double>::epsilon();

    // put the point back into the z-aligned aligned spheroid,
    // i.e. in referential R1
    UniformVectorType xs;
    for (unsigned int i=0; i<3; i++)
    {
      xs[i]=x[i];
    }
    
    xs[3] = 1;
    xs = m_InternalTransformInverse* xs;

    xsi[0] = std::sqrt (xs[0]*xs[0] + xs[1]*xs[1]);
    if ( (std::abs (xs[0]) < epsilon) && (std::abs (xs[1]) < epsilon) )
      xsi[1] = 0;
    else
    {
      if (xs[0] >= 0)
	xsi[1] = asin (xs[1] / xsi[0]);
      else
	xsi[1] = - asin (xs[1] / xsi[0]) + vnl_math::pi / 2.0;
    }
    xsi[2] = xs[2];

    // xsi has constrains and bounds, lets check them;
    if ( ( (xsi[0] < 0.0) ) ||
	 ( (xsi[1] < - vnl_math::pi / 2.0 ) || (xsi[1] > vnl_math::pi / 2.0) ) )
    {
      itkDebugMacro (<<"inconsistent Cylindrical Coordinate : "<<xsi);
    }

    // all singularity situations have been handled
    return xsi;
    
  }

  
  template <class TPixelType>
  typename CylindricalTransform<TPixelType>::VectorType CylindricalTransform<TPixelType>::ToCartesian(VectorType v, PointType p) const
  {
    
    // xsi has constrains and bounds, lets check them;
    if ( ( (p[0] < 0.0) ) ||
	 ( (p[1] < - vnl_math::pi / 2.0 ) || (p[1] > vnl_math::pi / 2.0) ) )
    {
      itkDebugMacro (<<"inconsistent Cylindrical Coordinate : "<<xsi);
    }
    
    MatrixType matrix; matrix = this->GetJacobianWithRespectToCoordinates(p);
    VectorType ret = matrix * v;

    return ret;
  }

  template <class TPixelType>
  typename CylindricalTransform<TPixelType>::VectorType CylindricalTransform<TPixelType>::ToCylindrical(VectorType v, PointType p) const
  {
    MatrixType matrix; matrix = this->GetJacobianWithRespectToCoordinates(p);
    MatrixType transpose; transpose = matrix.GetTranspose ();
    VectorType ret = transpose * v;
    
    return ret;
  }

  template <class TPixelType>
  void CylindricalTransform<TPixelType>::EvaluateScaleFactors (double xsi[3], double h[3]) const
  {
    // error underwhich we display a warning because we reached the singularity
    const double epsilon_singularity = 0.001;
    if ( xsi[0] <= epsilon_singularity )
    {
      itkWarningMacro (<<"singularity point : "<<xsi<<"\n"
		       <<"scale factors undefined --> set to null");
      h[0] = h[1] = h[2] = 0.0;
      return;
    }

    // norm of dX / dxsi[0]
    h[0] = 1.0 ;
    // norm of dX / dxsi[1]
    h[1] = 1.0 / xsi[0] ;
    // norm of dX / dxsi[2]
    h[2] = 1.0 ;
  }

  template <class TPixelType>
  void CylindricalTransform<TPixelType>::PrintSelf(std::ostream& os, Indent indent) const
  {
    this->Superclass::PrintSelf(os,indent);

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



  template <class TPixelType>
  double CylindricalTransform<TPixelType>::EstimateGeodesicLength1(PointType xi, VectorType dxi, unsigned int divisions) const
  {
    VectorType step = dxi / (double)(divisions);
    double length = 0.0;
    PointType point1 = xi, point2 = xi;
    PointType cart1, cart2;
    VectorType dL;
    
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
  double CylindricalTransform<TPixelType>::EstimateGeodesicLength2(PointType xi, VectorType dxi, unsigned int divisions) const
  {
    VectorType step = dxi / (double)(divisions);
    PointType pt = xi;
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
      if (pt[2] > (2 * vnl_math::pi))
        pt[2] -= 2 * vnl_math::pi;
      if (pt[2] < 0)
        pt[2] += 2 * vnl_math::pi;

    }

    return length;
  }


} // end namespace itk


#endif
