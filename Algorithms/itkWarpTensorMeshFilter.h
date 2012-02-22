/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkWarpTensorMeshFilter.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef __itkWarpTensorMeshFilter_h
#define __itkWarpTensorMeshFilter_h

#include <itkMeshToMeshFilter.h>
#include <itkMesh.h>
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkTensorLinearInterpolateImageFunction.h"
#include "itkWarpJacobianFilter.h"
#include "itkPoint.h"
#include "itkFixedArray.h"
#include "itkTensor.h"

#include <vector>

namespace itk
{

  template <class TMesh, class TDisplacementField>
    class ITK_EXPORT WarpTensorMeshFilter :
  public MeshToMeshFilter< TMesh, TMesh>
  {
  public:
    /** Standard class typedefs. */
    typedef WarpTensorMeshFilter Self;
    typedef TMesh MeshType;
    typedef MeshToMeshFilter< MeshType, MeshType > Superclass;

    typedef typename MeshType::Pointer    MeshPointer;
    typedef typename MeshType::PixelType  PixelType;

    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;
    
    itkNewMacro(Self);
    itkTypeMacro(WarpTensorMeshFilter, Superclass);

    typedef PixelType TensorType;
    typedef typename PixelType::ValueType ScalarType;
    typedef typename MeshType::PointType   PointType;

    enum ReorientationType { FS, PPD };

    /** Typedef to describe the region type. */
    typedef typename TDisplacementField::RegionType RegionType;

    /** Inherit some types from the superclass. */
    typedef TDisplacementField                           DisplacementFieldType;
    typedef typename DisplacementFieldType::Pointer      DisplacementFieldPointer;
    typedef typename DisplacementFieldType::ConstPointer DisplacementFieldConstPointer;
    typedef typename DisplacementFieldType::SpacingType  SpacingType;
    typedef typename DisplacementFieldType::IndexType    IndexType;
  
    /** Determine the image dimension. */
    itkStaticConstMacro(DisplacementFieldDimension, unsigned int,
			TDisplacementField::ImageDimension );

    itkStaticConstMacro(VectorDimension, unsigned int,
			TDisplacementField::PixelType::Dimension);

    typedef Matrix<double, DisplacementFieldDimension, VectorDimension> JacobianType;
    typedef Image<JacobianType, DisplacementFieldDimension>             JacobianImageType;
    typedef typename JacobianImageType::Pointer             JacobianPointer;
    typedef WarpJacobianFilter<TDisplacementField, JacobianImageType>  JacobianFilterType;

    /** Interpolator typedef support. */
    typedef VectorLinearInterpolateImageFunction<DisplacementFieldType,double>   DisplacementInterpolatorType;
    typedef typename DisplacementInterpolatorType::Pointer   DisplacementInterpolatorPointer;
  
    itkGetObjectMacro(Jacobian, JacobianImageType);

    /** Displacement field typedef support. */
    typedef typename DisplacementFieldType::PixelType DisplacementType;

    /** Interpolator typedef support. */
    typedef ScalarType CoordRepType;

    /** Type for representing the direction of the output image */
    typedef typename DisplacementFieldType::DirectionType     DirectionType;

    itkGetObjectMacro (DisplacementField, DisplacementFieldType);
    itkSetObjectMacro (DisplacementField, DisplacementFieldType);
    itkGetObjectMacro (InverseDisplacementField, DisplacementFieldType);
    itkSetObjectMacro (InverseDisplacementField, DisplacementFieldType);  

    /** Set the edge padding value */
    itkSetMacro( EdgePaddingValue, PixelType );

    /** Get the edge padding value */
    itkGetMacro( EdgePaddingValue, PixelType );
  
    /** Set the type of tensor reorientation */
    itkSetMacro( ReorientationStrategy, ReorientationType );

    /** Get the edge padding value */
    itkGetMacro( ReorientationStrategy, ReorientationType );

    void GenerateOutputInformation();

    void GenerateInputRequestedRegion();

  protected:
    WarpTensorMeshFilter();
    ~WarpTensorMeshFilter() {};
    void PrintSelf(std::ostream& os, Indent indent) const;

    /** WarpTensorMeshFilter is implemented as a filter.
     * As such, it needs to provide and implementation for 
     * GenerateData(). */
    void GenerateData();

  private:
    WarpTensorMeshFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    PixelType                  m_EdgePaddingValue;
    typename JacobianImageType::Pointer m_Jacobian;
    ReorientationType m_ReorientationStrategy;
    DisplacementInterpolatorPointer m_InverseDisplacementInterpolator;
    DisplacementFieldPointer m_DisplacementField;
    DisplacementFieldPointer m_InverseDisplacementField;
  
  
  };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWarpTensorMeshFilter.txx"
#endif

#endif
