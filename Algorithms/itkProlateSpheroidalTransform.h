/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkProlateSpheroidalTransform.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_ProlateSpheroidalTransform_h_
#define _itk_ProlateSpheroidalTransform_h_

#include "itkTransform.h"
#include "itkMatrix.h"
#include "itkPoint.h"
#include "itkImage.h"
#include "itkGradientImageFilter.h"
#include <itkLinearInterpolateImageFunction.h>
#include "itkCrossHelper.h"

#include <cmath>


namespace itk
{
  




  /**
     \class ProlateSpheroidalTransform (itk)
     \brief Transformation class from Cartesian coordinates to Prolate Spheroidal Coordinates.
     
     ProlateSpheroidalTransform is a transformation. Input parameters are 3 points defining a prolate
     spheroid: the center of the spheroid, the long axis point, and the short axis point.

     Use DefineEllipsoid() with the three points as arguments to define the spheroid (ellipsoid).   

     You can also use separately the methods SetCenter() SetLongAxisPoint() and SetShortAxisPoint().

     The ellipsoid is defined as :
     \f[

     \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2}  = 1

     \f]
     
     whith  \f$ a > b = c \f$. which means the ellipsoid is prolate and its long axis is the x-axis.

     There is an "internal" transformation InternalTransformInverse that transform a point (or vector)
     in space to the x-axis centered ellipsoid. We use usual definition of the prolate
     spheroidal coordinates :

     \f[
     
     \begin{array}{ccc}
     x = f_{1/2} \sinh({\xi}^1) \sin({\xi}^2) \cos({\xi}^3) \\
     y = f_{1/2} \sinh({\xi}^1) \sin({\xi}^2) \sin({\xi}^3) \\
     z = f_{1/2} \cosh({\xi}^1) \cos({\xi}^2)
     \end{array}
     
     \f]
     

     This class derives from an itk::Transform. Therefore beneficiates from all common methods for a
     transform, such as TransformPoint() or GetJacobian(). However, as this change of coordinates is not
     spacially homogeneous, you cannot use TransformVector() without indicating from which point this
     vector is defined.

     The coordinate change can be "Forward" - from Cartesian towards Prolate Spheroidal Coordinates - or
     "Backward" - from Prolate Spheroidal towards Cartesian coordinates. Control the direction of the transform
     with SetForward() methods. By definition, the inverse of a Forward transformation is the related Backward
     transformation. Therefore, Accessing GetInverse() has the same effect as using the same transform but changing
     direction to "Backward".     
     
     It can also be seen as an implicit function (EvaluateFunction()) using the formula:
     \f[

     f = \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} - 1

     \f]

     and its gradient (corresponding to first contravariant vector) using EvaluateGradient():
     \f[

     g = 

     \f]

     
     
     \ingroup Transform
     
     \sa ProlateSpheroidTransformTensorFilter
     \sa Transform
     
     \author Nicolas Toussaint, King's College London, nicolas.toussaint@kcl.ac.uk
  */
  
  template <class TPixelType=double>
  class ITK_EXPORT ProlateSpheroidalTransform : public Transform<TPixelType, 3, 3>
  {
  public:
    typedef ProlateSpheroidalTransform Self;
    typedef Transform<double, 3, 3>  Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;

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
    itkTypeMacro (ProlateSpheroidalTransform, Transform);

    /**
       Set the direction of the transformation. If Forward is set to 1, then the transformation is
       from Cartesian towards Prolate Spheroidal coordinates, and vice-versa.
    */
    itkSetClampMacro (Forward, unsigned int, 0, 1);
    itkGetConstMacro (Forward, unsigned int);
    itkBooleanMacro (Forward);

    /**
       Set the center point of the ellipsoid separately. Of course it has to be in Cartesian
       coordinates - as the center point has [0,0,0] as Prolate Spheroidal coordinates by definition.
    */
    virtual void SetCenter (PointType point)
    {
      m_Center = point;
      this->ComputeTransformation();
    }
    /**
       Set the long axis point of the ellipsoid separately. Of course it has to be in Cartesian
       coordinates - as the center point has [0,0,0] as Prolate Spheroidal coordinates by definition.
    */
    virtual void SetLongAxisPoint (PointType point)
    {
      m_LongAxisPoint = point;
      this->ComputeTransformation();
    }
    /**
       Set the short axis point of the ellipsoid separately. Of course it has to be in Cartesian
       coordinates.
    */
    virtual void SetShortAxisPoint (PointType point)
    {
      m_ShortAxisPoint = point;
      this->ComputeTransformation();
    }
    /**
       defining the prolate spheroid using three landmark points.
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
    virtual void DefineEllipsoid (PointType center,
			  PointType longaxispoint,
			  PointType shortaxispoint)
    {
      m_Center = center;
      m_LongAxisPoint = longaxispoint;
      m_ShortAxisPoint = shortaxispoint;
      this->ComputeTransformation();
    }
    
    /**  Method to transform a point. */
    virtual OutputPointType     TransformPoint(const InputPointType  & ) const;
    /**  Method to transform a vector. */
    virtual OutputVectorType    TransformVector(const InputVectorType &, const InputPointType &) const;
    /**  Method to transform a vnl_vector. */
    virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &, const InputPointType &) const;
    /**  Method to transform a CovariantVector. */
    virtual OutputCovariantVectorType TransformCovariantVector(const InputCovariantVectorType &, const InputPointType &) const;    
    /**
     * Set the transformation parameters and update internal transformation.
     * SetParameters gives the transform the option to set it's
     * parameters by keeping a reference to the parameters, or by
     * copying.  To force the transform to copy it's parameters call
     * SetParametersByValue.
     * 
     * In our case the parameters correspond to the cartesian coordinates of succesively
     * the center, the long axis and the long axis points of the ellipsoid.
     *
     * \sa SetParametersByValue
     * 
    */
    virtual void SetParameters( const ParametersType & );
    
    /**
     * Set the transformation parameters and update internal transformation. 
     * This method forces the transform to copy the parameters.  The
     * default implementation is to call SetParameters.  This call must
     * be overridden if the transform normally implements SetParameters
     * by keeping a reference to the parameters.
     * 
     * In our case the parameters correspond to the cartesian coordinates of succesively
     * the center, the long axis and the short axis points of the ellipsoid.
     * 
     * \sa SetParameters
     *
    */
    virtual void SetParametersByValue ( const ParametersType & p );

    /**
     * Get the Transformation Parameters.
     * In our case the parameters correspond to the cartesian coordinates of succesively
     * the center, the long axis and the long axis points of the ellipsoid.
     */
    virtual const ParametersType& GetParameters(void) const
    { 
      return this->m_Parameters; 
    }
    
    /** Set the fixed parameters and update internal transformation. */
    virtual void SetFixedParameters( const ParametersType &t ) 
    {
      if (t.GetSize())
	itkWarningMacro( << "There is no fixed parameters in a ProlateSpheroidalTransform" );
    }

    /** Get the Fixed Parameters. */
    virtual const ParametersType& GetFixedParameters(void) const
    {
      return this->m_FixedParameters;
    }
    
    /**
       Compute the Jacobian of the transformation
       *
       * This method computes the Jacobian matrix of the transformation
       * at a given input point. The rank of the Jacobian will also indicate 
       * if the transform is invertible at this point.
       *
       * The Jacobian is be expressed as a matrix of partial derivatives of the
       * output point components with respect to the parameters that defined
       * the transform:
       *
       * \f[
       *
       J=\left[ \begin{array}{cccc}
       \frac{\partial x_{1}}{\partial p_{1}} & 
       \frac{\partial x_{1}}{\partial p_{2}} & 
       \cdots  & \frac{\partial x_{1}}{\partial p_{m}}\	\
       \frac{\partial x_{2}}{\partial p_{1}} & 
       \frac{\partial x_{2}}{\partial p_{2}} & 
       \cdots  & \frac{\partial x_{2}}{\partial p_{m}}\	\
       \vdots  & \vdots  & \ddots  & \vdots \		\
       \frac{\partial x_{n}}{\partial p_{1}} & 
       \frac{\partial x_{n}}{\partial p_{2}} & 
       \cdots  & \frac{\partial x_{n}}{\partial p_{m}}
       \end{array}\right] 
       *
       * \f]
       */
    virtual vnl_matrix_fixed<TPixelType,3,3> GetJacobianWithRespectToCoordinates(const InputPointType  &) const;

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
       Get the inverse of this transform
    */
    virtual bool GetInverse( Self* inverse) const;
    
    /**
       Set the transform to identity. Basically remove (and unregister)
       the displacement field from the transform
    */
    virtual void SetIdentity(void);
    
    
    /**
       Indicates that this transform is linear. That is, given two
       points P and Q, and scalar coefficients a and b, then
       
       T( a*P + b*Q ) = a * T(P) + b * T(Q)
    */
    virtual bool IsLinear() const { return true; }
    
    /**
       from an input point \arg x (in Cartesian if Forward / in Prolate if !forward),
       
       Evaluate local basis at a defined prolate spheroidal coordinate \f$ \xi \f$ :

       \arg n : normal vector (radial) --> \f$ \xi^1 \f$.
       \arg a : right ascension vector --> \f$ \xi^2 \f$.
       \arg d : declination vector (circumferential) --> \f$ \xi^3 \f$.
       
       (n, a, d) is a local orthonormal prolate spheroid basis at position \f$\xi\f$. It is estimated
       using the contravariant vector formulae. For instance,
       \f$ n = g_1 =
       \left(
       \begin{array}{ccc}
       cosh({\xi}^1) sin({\xi}^2) cos({\xi}^3) \\
       cosh({\xi}^1) sin({\xi}^2) sin({\xi}^3) \\
       sinh({\xi}^1) cos({\xi}^2)
       \end{array}
       \right) 
       \f$

       \warning The contravariant basis is subject to a lot of sinus, cosinus, and hyperbolic function
       calculation, resulting in a possible numerical error, and inducing the non-orthogonality of the
       basis.       
    */
    virtual void EvaluateLocalBasis(InputPointType x,
			    VectorType &n,
			    VectorType &a,
			    VectorType &d) const;
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
    virtual PointType ToCartesian(PointType xsi) const;
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
    virtual PointType ToProlate(PointType x) const;

    virtual VectorType ToCartesian (VectorType v, PointType p) const;
    virtual VectorType ToProlate (VectorType v, PointType p) const;
    
    /// Recover the spheroid first eigenvalue
    itkGetMacro (Lambda1, double);
    /// Recover the spheroid second eigenvalue
    itkGetMacro (Lambda2, double);
    /// Recover the spheroid third eigenvalue
    itkGetMacro (Lambda3, double);
    /// Recover the spheroid eccentricity
    itkGetMacro (Eccentricity, double);    
    /// Recover the spheroid semi foci distance (a) 
    itkGetMacro (SemiFociDistance, double);
    /// Recover the spheroid first focus
    itkGetMacro (Focus1, PointType);
    /// Recover the spheroid second focus
    itkGetMacro (Focus2, PointType);
    /// Recover the spheroid center
    itkGetMacro (Center, PointType);
    /// Recover the first eigenvector
    itkGetMacro (Axis1, VectorType);
    /// Recover the second eigenvector
    itkGetMacro (Axis2, VectorType);
    /// Recover the third eigenvector
    itkGetMacro (Axis3, VectorType);
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
    virtual void EvaluateScaleFactors (double xsi[3], double h[3]) const;

    /**
       Evaluate ellipsoid equation. The implicit function that is evaluated here is :
       \f$
       f = \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} - 1
       \f$
       The argument is a Point in Cartesian coordinates
    */
    virtual double EvaluateFunction(PointType x) const;
    /**
       Evaluate ellipsoid equation. The implicit function that is evaluated here is :
       \f$
       f = \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} - 1
       \f$
       The argument is a Point in Cartesian coordinates
    */
    virtual double EvaluateFunction(double x[3]) const;
    /**
       Evaluate ellipsoid implicit function gradient. It is estimated in the centered ellipsoid
       coordinate system (point is transformed using InternalTransformInverse), with the coefficients
       of the ellipsoid. The argument is a Point in Cartesian coordinates \f$x\f$.
       \f[
       g_1 = 
       2 
       \left(
       \begin{array}{ccc}
       a^{-2} x(0) \\
       b^{-2} x(1) \\
       c^{-2} x(2)
       \end{array}
       \right)
       \f]

       It also corresponds to the first contravariant vector.

       \f[

       g_1 =
       \left(
       \begin{array}{ccc}
       cosh({\xi}^1) sin({\xi}^2) cos({\xi}^3) \\
       cosh({\xi}^1) sin({\xi}^2) sin({\xi}^3) \\
       sinh({\xi}^1) cos({\xi}^2)
       \end{array}
       \right) 

       \f]
    */
    virtual void EvaluateGradient(PointType x, VectorType &n) const;

    virtual double EstimateGeodesicLength1(PointType xi, VectorType dxi, unsigned int divisions) const;
    virtual double EstimateGeodesicLength2(PointType xi, VectorType dxi, unsigned int divisions) const;

  static double asinh(double value)
  { 
    double returned;
    
    if(value>0)
      returned = log(value + sqrt(value * value + 1));
    else
      returned = -log(-value + sqrt(value * value + 1));
    
    return(returned);
  }
  
  protected:

    ProlateSpheroidalTransform();
    ~ProlateSpheroidalTransform();
    
    virtual void PrintSelf(std::ostream& os, Indent indent) const;

    virtual void ComputeTransformation (void);
    
    PointType m_Center;
    PointType m_LongAxisPoint;
    PointType m_ShortAxisPoint;
    VectorType m_Axis1;
    VectorType m_Axis2;
    VectorType m_Axis3;

    double m_Coefficients[3];
    double m_Lambda1;
    double m_Lambda2;
    double m_Lambda3;

    PointType m_Focus1;
    PointType m_Focus2;
    double m_Eccentricity;
    double m_SemiFociDistance;
  
    UniformMatrixType m_InternalTransform;
    UniformMatrixType m_InternalTransformInverse;

    CrossHelperType m_Cross;
    unsigned int m_Forward;
    
  private :
    ProlateSpheroidalTransform(const Self&);
    void operator=(const Self&);


    
  };
  
}// end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkProlateSpheroidalTransform.txx"
#endif

#endif
