/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkTensorMeshCostFunction.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef __itkTensorMeshCostFunction_h
#define __itkTensorMeshCostFunction_h

#include "itkSingleValuedCostFunction.h"
#include "itkGaussianInterpolationTensorMeshFilter.h"
#include "itkWarpTensorMeshFilter.h"
#include "itkProlateSpheroidalTransformTensorMeshFilter.h"
#include <itkImage.h>
#include <itkBSplineKernelFunction.h>
#include <itkKaiserBesselKernelFunction.h>
#include <itkGaussianKernelFunction.h>

namespace itk
{
  
/** \class TensorMeshCostFunction
 *
 * \brief This class is a ITK CostFunction returning the difference between tensor fields
 *
 * This is a cost function between two tensor fields that are SPARSELY distributed in
 * space. This means the input of the cost function consists on two sets of tensors 
 * that are located in (possibly) different positions in space. The tensor sets do NOT
 * have to contain the same amount of tensors. 
 *
 * The mnost important method is SetVariables(). It takes as input the two tensor fields 
 * from the following arguments :
 * \param pts The first tensor field positions
 * \param tns The first tensor field tensors
 * \param refpts The second tensor field positions
 * \param reftns The second tensor field tensors
 *
 * The cost function is symmetric and it does not matter if you swap the 2 sets of tensors.
 *
 * **CAUTION** It is crucial that the amount of points matches the respective amount of
 * tensors.
 *
 * This cost function works in the Log-space. Thus the input tensors are "logged" prior
 * to processing. In order to be able to compare the two tensor fields, we need them to
 * be located at the same positions. Therefore a tri-linear anisotropic interpolation
 * is performed to the second tensor field onto the first tensor field positions. See 
 * InterpolateSparseTensors() for details.
 *
 * The interpolation uses a anisotropic Gaussian Kernel by default but you can modify
 * the kernel \f$ K_H\f$ through SetKernel(). 
 * 
 * The kernel estimate at position is given by:
 * 
 * \f{equation*}{
 * \hat{m}_{H}(D_{\xi}) = \exp \left( \frac{ \sum_{i=1}^N K_{H} (\xi - \xi_i) \log(D_{\xi_i})}{\sum_{i=1}^N
 * K_{H}(\xi -\xi_i)} \right)
 * \f}
 * 
 * The main purpose of this class is to return a estimation of the difference between
 * two input tensor fields given a set of input "parameters". The parameters are
 * the sigmas of the anisotropic tri-linear interpolation (Gaussian by default).
 *
 * If you work in Prolate Coordinate System, please consider using the
 * automatic coordinate switchers : SetChangeCoordinate() and don't forget to 
 * switch on SetPiWorkAround().
 *
 * \todo Give the possibility to the user to switch between angle difference and
 * Hosdorff distance between tensor fields.
 *
 *
 * \author Nicolas Toussaint, King's College London, nicolas.toussaint@kcl.ac.uk
 *
 * \ingroup TensorProcessing
 */
class ITK_EXPORT TensorMeshCostFunction : 
    public SingleValuedCostFunction 
{
public:
  /** Standard class typedefs. */
  typedef TensorMeshCostFunction              Self;
  typedef SingleValuedCostFunction            Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  typedef double                                                          ScalarType;
  typedef Tensor<ScalarType, 3>                                           TensorType;
  typedef DefaultStaticMeshTraits<TensorType, 3, 3, ScalarType,
    ScalarType, TensorType>                                               MeshTraits;
  typedef Mesh<TensorType, 3, MeshTraits>                                 TensorMeshType;
  typedef TensorMeshType::PointType                                       PointType;
  typedef TensorMeshType::PointIdentifier                                 PointIdentifier;
  
  typedef GaussianInterpolationTensorMeshFilter<TensorMeshType>           InterpolatorType;
  typedef Vector<ScalarType, 3>                                                VectorType;
  typedef VectorType                                                      DisplacementType;
  typedef Image<DisplacementType, 3>                                      DisplacementFieldType;
  typedef WarpTensorMeshFilter<TensorMeshType, DisplacementFieldType>     WarperType;
  typedef itk::ProlateSpheroidalTransformTensorMeshFilter<TensorMeshType> CoordinateSwitcherType;
  typedef CoordinateSwitcherType::TransformType                           TransformType;
  typedef TensorMeshType::PointsContainer                                 PointContainer;
  typedef TensorMeshType::PointDataContainer                              PixelContainer;

  typedef GaussianKernelFunction                                          GaussianKernelFunctionType;
  typedef BSplineKernelFunction<3>                                        BSplineKernelFunctionType;
  typedef KaiserBesselKernelFunction                                      KaiserBesselKernelType;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( TensorMeshCostFunction, SingleValuedCostFunction );
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  /**  MeasureType typedef.
   *  It defines a type used to return the cost function value. */
  typedef Superclass::MeasureType MeasureType;
  /** DerivativeType typedef.
   *  It defines a type used to return the cost function derivative.  */
  typedef Superclass::DerivativeType DerivativeType;
  /**  ParametersType typedef.
   *  It defines a position in the optimization search space. */
  typedef Superclass::ParametersType ParametersType;
  /**  ParametersType typedef.
   *  It defines a position in the optimization search space. */
  typedef Superclass::ParametersValueType ParametersValueType;
  /** Return the number of parameters required to compute 
   *  this cost function.
   *  This method MUST be overloaded by derived classes. */
  virtual unsigned int GetNumberOfParameters(void) const
  { return 3; }
  /** This method returns the value of the cost function corresponding
    * to the specified parameters.    */ 
  virtual MeasureType GetValue( const ParametersType & parameters ) const;
  /** This method returns the derivative of the cost function corresponding
    * to the specified parameters.   */ 
  virtual void GetDerivative( const ParametersType & parameters,
                              DerivativeType & derivative ) const;
  /** Get/Set the bounds of the cost function in the parameter space.
   *  If a value is asked with parameters outside those bounds, a maximal
   *  cost function value will be returned. */
  void GetBounds (ParametersValueType bds[2][100])
  {
    for (unsigned int i=0; i<2; i++)
      for (unsigned int j=0; j<100; j++)
	bds[i][j] = m_Bounds[i][j];
  }
  /** Get/Set the bounds of the cost function in the parameter space.
   *  If a value is asked with parameters outside those bounds, a maximal
   *  cost function value will be returned. */  
  void SetBounds (ParametersValueType bds[2][100])
  {
    for (unsigned int i=0; i<2; i++)
      for (unsigned int j=0; j<100; j++)
	m_Bounds[i][j] = bds[i][j];
  }

  itkGetMacro (UseProlateSpheroidalCoordinates, unsigned int);
  itkSetClampMacro (UseProlateSpheroidalCoordinates, unsigned int, 0, 1);
  itkBooleanMacro (UseProlateSpheroidalCoordinates);
  
  /** Set the variables of this cost function.
   *  \param pts is an array of positions where the measured data are.
   */
  void SetInputs( TensorMeshType::Pointer data,
		  TensorMeshType::Pointer domain,
		  TensorMeshType::Pointer reference);

  void SetDisplacementFields (DisplacementFieldType::Pointer forward,
			      DisplacementFieldType::Pointer backward)
  {
    m_ForwardDisplacementField  = forward;
    m_BackwardDisplacementField = backward;
    this->UpdatePipeline();
  }

  void SetTransform (TransformType::Pointer forward)
  {
    m_Transform = forward;
    this->UpdatePipeline();
  }

  bool Copy (const TensorMeshType::Pointer input, TensorMeshType::Pointer output) const;
  bool WriteMesh (TensorMeshType::Pointer mesh, const char* filename) const;

  enum KernelTypeIds
  {
    Gaussian = 0,
    BSpline,
    KaiserBessel
  };

  itkGetConstMacro (KernelType, unsigned int);
  itkSetClampMacro (KernelType, unsigned int, Gaussian, KaiserBessel);
  
protected:
  TensorMeshCostFunction();
  virtual ~TensorMeshCostFunction();
  
  /** This method root mean square distance (in the sense of Frobenius)
   *  of 2 distinct tensor fields
   *   Tensor arrays must have same size : no check is performed.    */ 
  virtual MeasureType EstimateMeanSquaredDistanceBetweenTensorFields( const TensorMeshType::Pointer tensors1,
								      const TensorMeshType::Pointer tensors2 ) const;
  /** This method outputs the mean angle (in degrees) between
   *  of 2 distinct tensor fields
   *   Tensor arrays must have same size : no check is performed.    */ 
  virtual MeasureType EstimateMeanAngleDifferenceBetweenTensorFields( const TensorMeshType::Pointer tensors1,
								      const TensorMeshType::Pointer tensors2 ) const;

  virtual void UpdatePipeline (void)
  {
    itkExceptionMacro (<<"should not enter this method. \n The metod should be overwritten ins subclasses");
  }
  
  TensorMeshType::Pointer m_RawData;
  TensorMeshType::Pointer m_InterpolatedData;
  TensorMeshType::Pointer m_Reference;  

  TransformType::Pointer m_Transform;
  DisplacementFieldType::Pointer m_ForwardDisplacementField;
  DisplacementFieldType::Pointer m_BackwardDisplacementField;
  
  bool m_UseProlateSpheroidalCoordinates;
  unsigned int m_KernelType;

  ParametersValueType m_WeightLimits[100];
  ParametersValueType m_Bounds[2][100];
  
private:
  TensorMeshCostFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
};

} // end namespace itk



#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorMeshCostFunction.cxx"
#endif

#endif
