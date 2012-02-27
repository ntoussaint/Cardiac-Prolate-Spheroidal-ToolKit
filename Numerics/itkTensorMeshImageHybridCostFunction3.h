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
#ifndef __itkTensorMeshImageHybridCostFunction3_h
#define __itkTensorMeshImageHybridCostFunction3_h

#include <itkSingleValuedCostFunction.h>

#include <itkBSplineKernelFunction.h>
#include <itkKaiserBesselKernelFunction.h>
#include <itkGaussianKernelFunction.h>

#include <itkTensorMeshToImageFilter.h>
#include <itkTensorImageToMeshFilter.h>
#include <itkLimitToAHAZoneImageFilter.h>
#include <itkWarpTensorMeshFilter.h>
#include <itkProlateSpheroidalTransformTensorMeshFilter.h>
#include <itkGaussianInterpolationTensorMeshFilter2.h>

namespace itk
{
  
class ITK_EXPORT TensorMeshImageHybridCostFunction3 : 
    public SingleValuedCostFunction 
{
public:
  /** Standard class typedefs. */
  typedef TensorMeshImageHybridCostFunction3  Self;
  typedef SingleValuedCostFunction            Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;
  itkTypeMacro( TensorMeshImageHybridCostFunction3, SingleValuedCostFunction );
  itkNewMacro(Self);

  typedef Superclass::MeasureType MeasureType;
  typedef Superclass::DerivativeType DerivativeType;
  typedef Superclass::ParametersType ParametersType;
  typedef Superclass::ParametersValueType ParametersValueType;

  typedef double                                 ScalarType;
  typedef ProlateSpheroidalTransform<ScalarType> TransformType;
  typedef Tensor<ScalarType, 3>                  TensorType;
  typedef Vector<ScalarType, 3>                       VectorType;
  typedef Image <VectorType, 3>                  DisplacementFieldType;
  typedef TensorImageToMeshFilter<TensorType, 3> ImageToMeshType;
  typedef TensorMeshToImageFilter<TensorType, 3> MeshToImageType;
  typedef ImageToMeshType::MeshType              MeshType;
  typedef ImageToMeshType::ImageType             ImageType;
  typedef MeshToImageType::DomainImageType       DomainImageType;
  typedef MeshType::PointType                    PointType;
  typedef MeshType::PointIdentifier              PointIdentifier;
  typedef MeshType::PointsContainer              PointContainer;
  typedef MeshType::PointDataContainer           PixelContainer;
  
  typedef WarpTensorMeshFilter<MeshType, DisplacementFieldType> WarperType;
  typedef ProlateSpheroidalTransformTensorMeshFilter<MeshType>  TransformerType;
  typedef GaussianInterpolationTensorMeshFilter2<MeshType>      Interpolator2Type;

  typedef GaussianKernelFunction                                GaussianKernelFunctionType;
  typedef BSplineKernelFunction<3>                              BSplineKernelFunctionType;
  typedef KaiserBesselKernelFunction                            KaiserBesselKernelType;
  
  enum KernelTypeIds
  {
    Gaussian = 0,
    BSpline,
    KaiserBessel
  };
  itkGetConstMacro (KernelType, unsigned int);
  itkSetClampMacro (KernelType, unsigned int, Gaussian, KaiserBessel);

  enum AttachTermTypeIds
  {
    MeanSquareDistance = 0,
    MeanSquareAngle
  };
  itkGetConstMacro (AttachTermType, unsigned int);
  itkSetClampMacro (AttachTermType, unsigned int, MeanSquareDistance, MeanSquareAngle);
  
  virtual unsigned int GetNumberOfParameters(void) const
  { return 3 * m_Interpolator->GetLimiter()->GetNumberOfAHAZones(); }
  
  virtual MeasureType GetValue      ( const ParametersType & parameters ) const;
  virtual void        GetDerivative ( const ParametersType & parameters, DerivativeType & derivative ) const
  {
    itkExceptionMacro (<<"NOT IMPLEMENTED");
  }
  
  
  virtual void GetOutputTerms (const ParametersType & parameters, double terms[2]) const;

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
  void SetInputs( MeshType::Pointer  data,
		  ImageType::Pointer reference)
  {
    m_InputData = data;
    m_InputReference = reference;
    this->UpdatePipeline();
  }
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
  
  /** Get/Set the lambda smoothness ratio to use as a trade-off between
   *  the data attach term and the regularity total variation term. */
  itkGetConstMacro(Lambda, double);
  /** Get/Set the lambda smoothness ratio to use as a trade-off between
   *  the data attach term and the regularity total variation term. */
  itkSetClampMacro(Lambda, double, 0, vcl_numeric_limits<double>::max());
  
  itkGetMacro (Beta, MeasureType);
  itkSetMacro (Beta, MeasureType);

  itkGetObjectMacro (Interpolator, Interpolator2Type);

protected:
  TensorMeshImageHybridCostFunction3();
  virtual ~TensorMeshImageHybridCostFunction3();

  MeshType::Pointer RecoverOutput (MeshType::Pointer mesh) const;
  
  virtual void GetValueInternal (const ParametersType & parameters, MeasureType terms[2]) const;  
  virtual void UpdatePipeline (void);

  /** This method mean square distance (in the sense of Frobenius)
   *  of 2 distinct tensor fields
   *   Tensor arrays must have same size : check is performed.    */ 
  virtual MeasureType EstimateMeanSquaredDistanceBetweenTensorFields( const MeshType::Pointer tensors1,
								      const MeshType::Pointer tensors2 ) const;
  /** This method  mean square angle (between main eigen vectors)
   *  of 2 distinct tensor fields
   *   Tensor arrays must have same size : check is performed.    */ 
  virtual MeasureType EstimateMeanSquaredAngleBetweenTensorFields( const MeshType::Pointer tensors1,
								      const MeshType::Pointer tensors2 ) const;
  virtual MeasureType EstimateTensorFieldTotalVariation (MeshType::Pointer mesh) const;
  
private:
  TensorMeshImageHybridCostFunction3(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  WarperType::Pointer        m_WarperIn;
  WarperType::Pointer        m_WarperOut;
  TransformerType::Pointer   m_TransformerIn;
  TransformerType::Pointer   m_TransformerOut;
  Interpolator2Type::Pointer m_Interpolator;
  
  MeshType::Pointer          m_InputData;
  ImageType::Pointer         m_InputReference;
  DomainImageType::Pointer   m_InputDomain;
  TransformType::Pointer     m_Transform;
  DisplacementFieldType::Pointer m_ForwardDisplacementField;
  DisplacementFieldType::Pointer m_BackwardDisplacementField;
  
  MeshType::Pointer    m_InternalData;
  MeshType::Pointer    m_InternalDomain;
  MeshType::Pointer    m_InternalReference;
  
  bool                m_UseProlateSpheroidalCoordinates;
  unsigned int        m_KernelType;
  ParametersValueType m_Bounds[2][100];
  MeasureType         m_Lambda;
  MeasureType         m_Beta;

  unsigned int m_AttachTermType;

  void Copy (const MeshType::Pointer input, MeshType::Pointer output) const;
  ImageType::Pointer CopyMeshToImage (MeshType::Pointer mesh, DomainImageType::Pointer domain) const;
  DomainImageType::Pointer CreateDomain (ImageType::Pointer image) const;

};

} // end namespace itk



#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorMeshImageHybridCostFunction3.cxx"
#endif

#endif
