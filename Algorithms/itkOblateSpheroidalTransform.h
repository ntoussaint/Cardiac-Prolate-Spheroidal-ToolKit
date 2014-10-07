/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkOblateSpheroidalTransform.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_OblateSpheroidalTransform_h_
#define _itk_OblateSpheroidalTransform_h_

#include "itkProlateSpheroidalTransform.h"

namespace itk
{
  /**
     \class OblateSpheroidalTransform (itk)
     \brief Transformation class from Cartesian coordinates to Oblate Spheroidal Coordinates.
     
     
     \ingroup Transform
     
     \sa ProlateSpheroidTransformTensorFilter
     \sa Transform
     
     \author Nicolas Toussaint, King's College London, nicolas.toussaint@kcl.ac.uk
  */
  
  template <class TPixelType=double>
  class ITK_EXPORT OblateSpheroidalTransform : public ProlateSpheroidalTransform<TPixelType>
  {
  public:
    typedef OblateSpheroidalTransform Self;
    typedef ProlateSpheroidalTransform<double>  Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;
    typedef SmartPointer<Superclass> SuperClassPointer;

    /** Type of the scalar representing coordinate and vector elements. */
    typedef TPixelType ScalarType;
    
    /** Dimension of the domain space. */
    itkStaticConstMacro(InputSpaceDimension, unsigned int,  3);
    itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(ParametersDimension, unsigned int,  9);
  
    /** Type of the input parameters. */
    typedef Superclass::ParametersType ParametersType;
    
    /** Type of the Jacobian matrix. */
    typedef Superclass::JacobianType JacobianType;
    
    /** Standard vector type for this class. */
    typedef Superclass::InputVectorType InputVectorType;
    typedef Superclass::OutputVectorType OutputVectorType;
    
    /** Standard covariant vector type for this class */
    typedef Superclass::InputCovariantVectorType InputCovariantVectorType;
    typedef Superclass::OutputCovariantVectorType OutputCovariantVectorType;
  
    /** Standard vnl_vector type for this class. */
    typedef Superclass::InputVnlVectorType InputVnlVectorType;
    typedef Superclass::OutputVnlVectorType OutputVnlVectorType;
    
    /** Standard coordinate point type for this class */
    typedef Superclass::InputPointType InputPointType;
    typedef Superclass::OutputPointType OutputPointType;

    /** Base inverse transform type. This type should not be changed to the
     * concrete inverse transform type or inheritance would be lost. */
    typedef Superclass::InverseTransformBaseType InverseTransformBaseType;
    typedef InverseTransformBaseType::Pointer    InverseTransformBasePointer;

    /// typedefs definitions
    typedef Matrix<ScalarType, 4, 4> UniformMatrixType;
    /// typedefs definitions
    typedef Matrix<ScalarType, 3, 3> MatrixType;
    /// typedefs definitions
    typedef Vector<ScalarType, 4> UniformVectorType;
    /// typedefs definitions
    typedef Point<ScalarType, 3> PointType;
    /// typedefs definitions
    typedef Vector<ScalarType, 3> VectorType;
    /// typedefs definitions
    typedef itk::CrossHelper<VectorType> CrossHelperType;
  
    /** generic constructors */
    itkNewMacro (Self);
    itkTypeMacro (OblateSpheroidalTransform, Transform);
  
    /**
       defining the oblate spheroid using four landmark points.
       The center of the spheroid, the long axis point and the short axis point.
       From there a ellipsoid is defined using the expression
       \f$
       \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2}  = 1
       \f$
       
       The spheroid is centered at origin [0,0,0] and its long axis to the x-axis, with \f$a > b = c\f$.
       The transformation from a point in space to the centered spheroid is simply a translation
       (using the center point) and a rotation (orthonormal matrix), stored in m_InternalTransformInverse.

       Some parameters are derived from this definition, such as :
       \arg The eigenvectors of the spheroid, normalized : \f$ ev_1 = longaxispoint - center \f$; \f$ ev_2 = shortaxispoint - center \f$;
       and \f$ ev_3 = ev_1 \times ev_2 \f$.
       \arg The eigenvalues of the spheroid taken as norms of the eigenvectors (obviously before normalization).
       \arg The coefficients of the ellipsoid: \f$ a = \|ev_1\| \f$, and \f$ b = c = \|ev_2\| \f$.
       \arg The semi-foci distance \f$f_{1/2}\f$, half of distance between the 2 foci, is also the linear eccentricity \f$l = \sqrt{a^2 - b^2}}\f$.
       \arg The eccentricity : \f$ e = \frac{\sqrt{a^2 - b^2}}{a} = l / a\f$.
    */
    void DefineOblateSpheroid (PointType center,
			       PointType longaxispoint,
			       PointType shortaxispoint)
    {
      this->DefineEllipsoid (center, longaxispoint, shortaxispoint );
    }
  
    /**
       Get the jacobian of the transformation with respect to the coordinates at a specific point
       
       It can also be seen in our case as a local contravariant basis, 
    */
    vnl_matrix_fixed<TPixelType,3,3> GetJacobianWithRespectToCoordinates(const InputPointType  &) const;
    /** Compute the Jacobian Matrix of the transformation at one point,
     *  allowing for thread-safety. */
    virtual void ComputeJacobianWithRespectToParameters( const InputPointType  & p, JacobianType & jacobian) const
    {  
      itkExceptionMacro( "ComputeJacobianWithRespectToPosition not yet implemented "
                         "for " << this->GetNameOfClass() );
    }
  
    virtual void ComputeJacobianWithRespectToPosition(const InputPointType &,
                                                      JacobianType &) const
    {  
      itkExceptionMacro( "ComputeJacobianWithRespectToPosition not yet implemented "
                         "for " << this->GetNameOfClass() );
    }
    
    /**
       Transforms a point between prolate spheroidal coordinates towards Cartesian coordinates.
       \f[
       
       \begin{array}{ccc}
       x(0) = f_{1/2} \sinh({\xi}^1) \sin({\xi}^2) \cos({\xi}^3) \\
       x(1) = f_{1/2} \sinh({\xi}^1) \sin({\xi}^2) \sin({\xi}^3) \\
       x(2) = f_{1/2} \cosh({\xi}^1) \cos({\xi}^2)
       \end{array}
       
       \f]

       and then the transformation InternalTransform iis used to go from zero-centered ellipsoid
       to the true point in space.
    */
    PointType ToCartesian(PointType xsi) const;
    /**
       Not as easy as it seems...

       By using the trignometric identity \f$\cos^2(x) + \sin^2(x) = 1\f$ and the hyperbolic identity
       \f$\cosh^2(x) - \sinh^2(x) = 1\f$, we find \f$ x(0)^2 + x(1)^2 = f_{1/2}^2 \sinh^2(\xi^1) \sin^2(\xi^2) \f$.

       If we note \f$A = f_{1/2}^2 \f$; \f$ B = x(0)^2 + x(1)^2 \f$; \f$ C = x(2)^2 \f$; and \f$ \alpha = \sin^2(\xi^2) \f$.
       Then the hyperbolic identity gives us a polynome in \f$\alpha\f$ :
       \f[
       A \alpha^2 + (-A + B + C) \alpha - B = 0
       \f]

       Since only one of the root \f$\Gamma_1\f$ is possible, we find \f$\xi^2\f$ (\f$\xi^2\f$ is positive).
       then from the trigonometric identity, we have \f$ \sin^2(\xi^1) = B / (A sin^2(\xi^2))\f$, from which we
       find \f$\xi^1\f$ since also positive.

       \f$\xi^3\f$ is simply derived from the projection of the cartesian point onto the plane \f$(x(1), x(2))\f$, i.e. \f$\xi^3 = \arccos \left( \frac{x(1)}{\|x\|} \right) \f$.
    */
    PointType ToProlate(PointType x) const;

    /**
       evaluate scale factors at prolate position \f$xi\f$:
       \f$ 
       \left(
       \begin{array}{cc}
       \frac{1}{h_1} = \frac{1}{h_2} = f_{1/2} \sqrt( \left( \sinh^2(\xi_1) + \sin^2(\xi_2) \right) 
       \frac{1}{h_3} = f_{1/2} \sinh(\xi_1) \sin(\xi_2)
       \end{array}
       \right) 
       \f$
    */
    void EvaluateScaleFactors (double xsi[3], double h[3]) const;

  protected:

    OblateSpheroidalTransform();
    ~OblateSpheroidalTransform();
    
    void ComputeTransformation (void);

  private :
    OblateSpheroidalTransform(const Self&);
    void operator=(const Self&);


    
  };
  
}// end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOblateSpheroidalTransform.txx"
#endif

#endif
