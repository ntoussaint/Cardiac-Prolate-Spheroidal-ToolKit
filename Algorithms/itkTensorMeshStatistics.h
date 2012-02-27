/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkTensorMeshStatistics.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_TensorMeshStatistics_h_
#define _itk_TensorMeshStatistics_h_

#include <itkImageToImageFilter.h>
#include <itkMesh.h>
#include <itkPoint.h>

#include <itkTensor.h>
#include <itkTensorImageToMeshFilter.h>
#include <itkTensorMeshToImageFilter.h>
#include <itkLimitToAHAZoneImageFilter.h>
#include "itkProlateSpheroidalTransformTensorMeshFilter.h"
#include "itkWarpTensorMeshFilter.h"

#include <vector>

namespace itk
{

  /**
     \class TensorMeshStatistics
     \brief 
     Author: Nicolas Toussaint. Copyright KCL 2010.
  */
  template < class TPrecision, unsigned int TDimension>
  class ITK_EXPORT TensorMeshStatistics :
    public ImageToImageFilter< Image < Tensor<TPrecision,TDimension>, TDimension >, Image < TPrecision,TDimension > >
  {
    
  public:

    typedef TensorMeshStatistics Self;
    typedef ImageToImageFilter< Image < Tensor<TPrecision,TDimension>, TDimension >, Image < TPrecision,TDimension > > Superclass;
    
    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;
    
    itkNewMacro(Self);
    itkTypeMacro(TensorMeshStatistics, ImageToImageFilter);

    typedef TPrecision                             ScalarType;
    typedef Tensor<ScalarType,TDimension>          TensorType;
    typedef ProlateSpheroidalTransform<ScalarType> TransformType;
    typedef Image<TensorType, TDimension>          InputImageType;
    typedef Image<TPrecision, TDimension>          OutputImageType;
    typedef Vector<ScalarType, TDimension>              DisplacementType;
    typedef Image<DisplacementType, TDimension>    DisplacementFieldType;

    typedef LimitToAHAZoneImageFilter<InputImageType>             AHALimiterType;
    typedef TensorImageToMeshFilter<TensorType, TDimension>       ImageToMeshFilterType;
    typedef TensorMeshToImageFilter<TensorType, TDimension>       MeshToImageFilterType;
    typedef typename ImageToMeshFilterType::MeshType              MeshType;
    typedef typename ImageToMeshFilterType::ImageType             ImageType;
    typedef typename MeshToImageFilterType::DomainImageType       DomainImageType;
    typedef ProlateSpheroidalTransformTensorMeshFilter<MeshType>  TransformerType;    
    typedef WarpTensorMeshFilter<MeshType, DisplacementFieldType> WarperType;

    typedef Tensor<ScalarType, TensorType::DegreesOfFreedom>      CovarianceMatrixType;
    typedef Vector<ScalarType, TensorType::DegreesOfFreedom>      TensorVectorType;

    typedef typename TensorType::VectorType VectorType;
    typedef typename MeshType::PointType PointType;
    
    enum CoordinatesTypesIds
    {
      ProlateSpheroidal = 0,
      Cartesian
    };
    itkGetMacro (CoordinatesType, unsigned int);
    itkSetClampMacro (CoordinatesType, unsigned int, ProlateSpheroidal, Cartesian);

    enum StatisticsOutputTypeIds
    {
      CovarianceMatrixNorm = 0,
      CovarianceMatrixFA,
      VarianceAroundEV1
    };
    
    itkGetMacro (StatisticsOutputType, unsigned int);
    itkSetClampMacro (StatisticsOutputType, unsigned int, CovarianceMatrixNorm, VarianceAroundEV1);
    
    itkGetObjectMacro (Transform, TransformType);
    itkSetObjectMacro (Transform, TransformType);
    itkGetObjectMacro (DisplacementField, DisplacementFieldType);
    itkSetObjectMacro (DisplacementField, DisplacementFieldType);
    itkGetObjectMacro (InverseDisplacementField, DisplacementFieldType);
    itkSetObjectMacro (InverseDisplacementField, DisplacementFieldType);

    itkGetObjectMacro (MeshOutput, MeshType);

    itkGetObjectMacro (Limiter, AHALimiterType);
    
  protected:

    TensorMeshStatistics();
    ~TensorMeshStatistics(){};
    
    virtual void GenerateData (void);
    virtual void GenerateOutputInformation(void);
    
    CovarianceMatrixType ComputeCovarianceMatrix (typename MeshType::Pointer mesh);
    TensorType ComputeDispersionTensor (typename MeshType::Pointer mesh);
    TensorType ComputeDispersionTensor2 (typename MeshType::Pointer mesh);
    TensorType ComputeGradientTensor (typename MeshType::Pointer mesh);
    ScalarType EvaluateScalar(CovarianceMatrixType sigma);

    typename DomainImageType::Pointer CreateDomain (typename MeshType::Pointer mesh);
    typename DomainImageType::Pointer CreateDomain (typename ImageType::ConstPointer input);
    
    TensorVectorType Tensor2Vec(const TensorType &tensor)
    {
      TensorVectorType vec;
      for( unsigned int i=0; i<TensorType::DegreesOfFreedom; i++)
	vec[i] = tensor.GetNthComponentAsVector(i);
      return vec;
    }

    void WriteTensorMesh  (typename MeshType::Pointer mesh,   const char* filename) const;
    void WriteTensorImage (typename ImageType::Pointer image, const char* filename) const;
    
    typename TransformType::Pointer          m_Transform;
    typename DisplacementFieldType::Pointer  m_DisplacementField;
    typename DisplacementFieldType::Pointer  m_InverseDisplacementField;
    typename AHALimiterType::Pointer         m_Limiter;
    typename ImageToMeshFilterType::Pointer  m_ImageToMeshFilter;
    typename WarperType::Pointer             m_Warper;
    typename WarperType::Pointer             m_InverseWarper;
    typename TransformerType::Pointer        m_Transformer;
    typename TransformerType::Pointer        m_InverseTransformer;
    typename MeshType::Pointer               m_AuxMesh;
    typename TransformerType::Pointer        m_AuxTransformer;
    typename WarperType::Pointer             m_AuxWarper;
    typename MeshType::Pointer               m_MeshOutput;

    unsigned int m_StatisticsOutputType;
    unsigned int m_CoordinatesType;
    
  private:
    TensorMeshStatistics(const Self&); // purposely not implemented
    void operator=(const Self&); // purposely not implemented

  };

} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorMeshStatistics.txx"
#endif

#endif
