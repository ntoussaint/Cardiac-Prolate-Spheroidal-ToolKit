/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkProlateSpheroidalTransformTensorMeshFilter.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_ProlateSpheroidalTransformTensorMeshFilter_h_
#define _itk_ProlateSpheroidalTransformTensorMeshFilter_h_


#include "itkMeshToMeshFilter.h"
#include "itkMesh.h"

#include "itkTensor.h"
#include "itkProlateSpheroidalTransform.h"

namespace itk
{
  /*! \class ProlateSpheroidalTransformTensorMeshFilter
    \ingroup 
  */
  
  template <class TMesh = Mesh < Tensor<double,3>, 3> > 
    class ITK_EXPORT ProlateSpheroidalTransformTensorMeshFilter :
    public MeshToMeshFilter< TMesh, TMesh>
    {
      
  public:

    typedef ProlateSpheroidalTransformTensorMeshFilter Self;
    typedef TMesh MeshType;
    typedef MeshToMeshFilter< MeshType, MeshType > Superclass;
  
    typedef typename MeshType::Pointer    MeshPointer;
    typedef typename MeshType::PixelType  PixelType;

    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;
    
    itkNewMacro(Self);
    itkTypeMacro(ProlateSpheroidalTransformTensorMeshFilter, Superclass);
 
    /** Standard class typedefs */

    /** Image typedefs */
    typedef PixelType TensorType;
    typedef typename PixelType::ValueType ScalarType;
    typedef typename MeshType::PointType   PointType;
    typedef Vector <ScalarType, 3>   VectorType;
    typedef ScalarType  RealType;
    typedef Matrix<RealType, 3, 3> MatrixType;
    typedef ProlateSpheroidalTransform<ScalarType> TransformType;
    typedef typename TransformType::Pointer TransformPointerType;

    void SetTransform( TransformType* transform );
    itkGetConstObjectMacro( Transform, TransformType );

    unsigned int GetForward (void)
    { return m_Transform->GetForward(); }
    void SetForward (unsigned int f)
    { m_Transform->SetForward (f); }
    itkBooleanMacro (Forward);

    void GenerateOutputInformation();
    
  protected:

    ProlateSpheroidalTransformTensorMeshFilter();
    virtual ~ProlateSpheroidalTransformTensorMeshFilter(){}
    
    void GenerateData ();
    void PrintSelf (std::ostream&, Indent) const;
    
  private:
    ProlateSpheroidalTransformTensorMeshFilter (const Self&);
    void operator=(const Self&);

    TransformPointerType m_Transform;
    
  };
  


} // end of namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkProlateSpheroidalTransformTensorMeshFilter.txx"
#endif

#endif
