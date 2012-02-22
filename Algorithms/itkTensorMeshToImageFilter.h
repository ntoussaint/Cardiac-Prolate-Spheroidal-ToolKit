/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkTensorMeshToImageFilter.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_TensorMeshToImageFilter_h_
#define _itk_TensorMeshToImageFilter_h_

#include <itkImageSource.h>
#include <itkMesh.h>
#include <itkTensor.h>

namespace itk
{  
  template <class TPixel = Tensor<double,3>, unsigned int TDimension = 3>
    class ITK_EXPORT TensorMeshToImageFilter :
    public ImageSource< Image<TPixel, TDimension> >
    {

    public:
    
    typedef TensorMeshToImageFilter Self;
    typedef ImageSource<Image<TPixel, TDimension> > Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    typedef TPixel                                        PixelType;
    typedef Image<TPixel, TDimension>                     ImageType;
    typedef typename PixelType::ValueType                 ScalarType;
    typedef DefaultStaticMeshTraits<PixelType, TDimension, TDimension, ScalarType, ScalarType, PixelType> MeshTraits;
    typedef Mesh<PixelType, TDimension, MeshTraits>       MeshType;
    typedef Image<ScalarType, TDimension>                 DomainImageType;
  
    itkNewMacro(Self);
    itkTypeMacro(TensorMeshToImageFilter, ImageSource);

    itkGetConstObjectMacro (Domain, DomainImageType);
    itkSetObjectMacro (Domain, DomainImageType);
    
    virtual void SetInput( unsigned int, const MeshType* mesh);
    virtual void SetInput( const MeshType* mesh)
    {
      this->SetInput (0, mesh);
    }
    const MeshType* GetInput (void)
    {
      return dynamic_cast< MeshType *>( this->ProcessObject::GetInput(0) );
    }
  
    protected:
    TensorMeshToImageFilter();
    ~TensorMeshToImageFilter(){};

    virtual void GenerateData(void);
    virtual void GenerateOutputInformation();

    virtual void PrintSelf(std::ostream& os, Indent indent) const
    {
      Superclass::PrintSelf(os,indent);
    }

    private:
    
    typename MeshType::Pointer m_Input;
    typename DomainImageType::Pointer m_Domain;
  
    };

} // end of namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorMeshToImageFilter.txx"
#endif

#endif
