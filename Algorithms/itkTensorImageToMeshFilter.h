/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkTensorImageToMeshFilter.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_TensorImageToMeshFilter_h_
#define _itk_TensorImageToMeshFilter_h_

#include <itkMeshSource.h>
#include <itkMesh.h>
#include <itkTensor.h>

namespace itk
{
  
  template <class TPixel = Tensor<double,3>, unsigned int TDimension = 3>
    class ITK_EXPORT TensorImageToMeshFilter :
    public MeshSource< Mesh <TPixel, TDimension, DefaultStaticMeshTraits <TPixel, TDimension, TDimension, typename TPixel::ValueType, typename TPixel::ValueType, TPixel> > >
    {

    public:
    
    typedef TensorImageToMeshFilter Self;
    typedef MeshSource< Mesh<TPixel, TDimension>  > Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    typedef TPixel                                        PixelType;
    typedef Image<TPixel, TDimension>                     ImageType;
    typedef typename PixelType::ValueType                 ScalarType;
    typedef DefaultStaticMeshTraits<PixelType, TDimension, TDimension, ScalarType, ScalarType, PixelType> MeshTraits;
    typedef Mesh<PixelType, TDimension, MeshTraits>       MeshType;
    typedef typename MeshType::PointType                  PointType;

    itkNewMacro(Self);
    itkTypeMacro(TensorImageToMeshFilter, MeshSource);
    itkStaticConstMacro( Dimension, unsigned int, TDimension);

    virtual void SetInput( unsigned int, const ImageType* image);
    virtual void SetInput( const ImageType* image)
    {
      this->SetInput (0, image);
    }
    const ImageType* GetInput (void)
    {
      return dynamic_cast< ImageType *>( this->ProcessObject::GetInput(0) );
    }
    
    protected:
    TensorImageToMeshFilter();
    ~TensorImageToMeshFilter(){};

    virtual void GenerateData(void);
    virtual void GenerateOutputInformation();

    private:

    typename ImageType::Pointer m_Input;

    };


} // end of namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorImageToMeshFilter.txx"
#endif

#endif
