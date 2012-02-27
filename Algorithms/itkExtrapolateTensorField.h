/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkExtrapolateTensorField.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_ExtrapolateTensorField_h_
#define _itk_ExtrapolateTensorField_h_

#include <itkMeshToMeshFilter.h>

#include <itkTensorImageToMeshFilter.h>
#include <itkTensorMeshToImageFilter.h>
#include <itkWarpTensorMeshFilter.h>

#include <itkProlateSpheroidalTransform.h>
#include <itkProlateSpheroidalTransformTensorMeshFilter.h>
#include <itkGaussianInterpolationTensorMeshFilter2.h>

#include <itkBSplineKernelFunction.h>
#include <itkKaiserBesselKernelFunction.h>
#include <itkGaussianKernelFunction.h>

namespace itk
{

  /**
     \class ExtrapolateTensorField
     \brief 
     Author: Nicolas Toussaint. Copyright KCL 2010.
  */
  template < class TPrecision, unsigned int TDimension>
  class ITK_EXPORT ExtrapolateTensorField :
    public MeshToMeshFilter< Mesh <Tensor<TPrecision,TDimension>, TDimension, DefaultStaticMeshTraits <Tensor<TPrecision,TDimension>, TDimension, TDimension, typename Tensor<TPrecision,TDimension>::ValueType, typename Tensor<TPrecision,TDimension>::ValueType, Tensor<TPrecision,TDimension> > >, Mesh <Tensor<TPrecision,TDimension>, TDimension, DefaultStaticMeshTraits <Tensor<TPrecision,TDimension>, TDimension, TDimension, typename Tensor<TPrecision,TDimension>::ValueType, typename Tensor<TPrecision,TDimension>::ValueType, Tensor<TPrecision,TDimension> > > >
  {
    
  public:

    typedef ExtrapolateTensorField Self;
    typedef MeshToMeshFilter< Mesh <Tensor<TPrecision,TDimension>, TDimension, DefaultStaticMeshTraits <Tensor<TPrecision,TDimension>, TDimension, TDimension, typename Tensor<TPrecision,TDimension>::ValueType, typename Tensor<TPrecision,TDimension>::ValueType, Tensor<TPrecision,TDimension> > >, Mesh <Tensor<TPrecision,TDimension>, TDimension, DefaultStaticMeshTraits <Tensor<TPrecision,TDimension>, TDimension, TDimension, typename Tensor<TPrecision,TDimension>::ValueType, typename Tensor<TPrecision,TDimension>::ValueType, Tensor<TPrecision,TDimension> > > > Superclass;
    
    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;
    
    itkNewMacro(Self);
    itkTypeMacro(ExtrapolateTensorField, ImageToImageFilter);

    typedef TPrecision                             ScalarType;
    typedef Tensor<ScalarType,TDimension>          TensorType;
    typedef ProlateSpheroidalTransform<ScalarType> TransformType;
    typedef Image<TPrecision, TDimension>          ImageType;
    typedef Vector<ScalarType, TDimension>              DisplacementType;
    typedef Image<DisplacementType, TDimension>    DisplacementFieldType;

    typedef TensorImageToMeshFilter<TensorType, TDimension>       ImageToMeshFilterType;
    typedef typename ImageToMeshFilterType::MeshType              MeshType;
    typedef ProlateSpheroidalTransformTensorMeshFilter<MeshType>  TransformerType;    
    typedef WarpTensorMeshFilter<MeshType, DisplacementFieldType> WarperType;
    typedef itk::GaussianInterpolationTensorMeshFilter2<MeshType>  InterpolatorType;
    
    itkGetObjectMacro (Transform, TransformType);
    itkSetObjectMacro (Transform, TransformType);
    itkGetObjectMacro (DisplacementField, DisplacementFieldType);
    itkSetObjectMacro (DisplacementField, DisplacementFieldType);
    itkGetObjectMacro (InverseDisplacementField, DisplacementFieldType);
    itkSetObjectMacro (InverseDisplacementField, DisplacementFieldType);
    itkGetObjectMacro (Domain, ImageType);
    itkSetObjectMacro (Domain, ImageType);
    
    itkGetObjectMacro (Interpolator, InterpolatorType);
    
    void SetAlpha (double* alpha)
    { m_Interpolator->SetAlpha (alpha); }
    void SetKernel (typename InterpolatorType::KernelFunctionType* kernel)
    { m_Interpolator->SetKernel (kernel); }
    void SetUseProlateCoordinates (unsigned int val)
    { m_UseProlateCoordinates = val; m_Interpolator->SetUsePiWorkAround (val); }
    itkBooleanMacro (UseProlateCoordinates);
    
  protected:

    ExtrapolateTensorField();
    ~ExtrapolateTensorField(){};
    
    virtual void GenerateData (void);
    virtual void GenerateOutputInformation(void);

    typename MeshType::Pointer DomainToMesh (typename ImageType::Pointer image);
    
    typename ImageType::Pointer              m_Domain;
    typename TransformType::Pointer          m_Transform;
    typename DisplacementFieldType::Pointer  m_DisplacementField;
    typename DisplacementFieldType::Pointer  m_InverseDisplacementField;
    typename WarperType::Pointer             m_Warper;
    typename WarperType::Pointer             m_DomainWarper;
    typename WarperType::Pointer             m_InverseWarper;
    typename TransformerType::Pointer        m_Transformer;
    typename TransformerType::Pointer        m_DomainTransformer;
    typename TransformerType::Pointer        m_InverseTransformer;
    typename InterpolatorType::Pointer       m_Interpolator;

    unsigned int m_UseProlateCoordinates;
    
  private:
    ExtrapolateTensorField(const Self&); // purposely not implemented
    void operator=(const Self&); // purposely not implemented

  };

} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkExtrapolateTensorField.txx"
#endif

#endif
