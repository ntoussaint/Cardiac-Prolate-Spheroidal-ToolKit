/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkProlateSpheroidalGradientTensorMeshFilter.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_ProlateSpheroidalGradientTensorMeshFilter_h_
#define _itk_ProlateSpheroidalGradientTensorMeshFilter_h_

#include <itkMeshToMeshFilter.h>
#include <itkMesh.h>
#include <itkTensor.h>
#include <itkPoint.h>
#include <itkProlateSpheroidalTransform.h>

#include <vector>

namespace itk
{

  /**
     \class ProlateSpheroidalGradientTensorMeshFilter
     \brief This class is a kernel based interpolation for tensors embedded in a mesh grid

     ProlateSpheroidalGradientTensorMeshFilter is a mesh to mesh filter that takes two inputs.
     Use SetInput(0, ...) to give the data to be used for interpolation. Use SetInput(1, ...)
     to give the mesh points where the evaluation is performed. The output mesh has the same
     points and structure as the second input (sampling points). At each point position of the
     output mesh, 

     In case of Non-Cartesian coordinates spatial positions (spherical, cylindrical, or prolate
     spheroidal), you should use UsePiWorkAroundOn(). This method will trigger a check during
     the spatial distance estimation, and module around -pi and pi for the second and third component,
     as the first one is supposed to be a radial distance.
     
     \sa GaussianKernelFunction MeshToMeshFilter
     Author: Nicolas Toussaint. Copyright KCL 2012.
  */
  
  template <class TMesh = Mesh < Tensor<double,3>, 3> >
  class ITK_EXPORT ProlateSpheroidalGradientTensorMeshFilter :
  public MeshToMeshFilter< TMesh, TMesh>
  {
    
  public:

  typedef ProlateSpheroidalGradientTensorMeshFilter Self;
    typedef TMesh MeshType;
    typedef MeshToMeshFilter< MeshType, MeshType > Superclass;

    typedef typename MeshType::Pointer    MeshPointer;
    typedef typename MeshType::PixelType  PixelType;

    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;
    
    itkNewMacro(Self);
    itkTypeMacro(ProlateSpheroidalGradientTensorMeshFilter, Superclass);
  
  /** Image typedefs */
    typedef PixelType                     TensorType;
    typedef typename PixelType::ValueType ScalarType;
    typedef typename MeshType::PointType  PointType;
    typedef Vector <ScalarType, 3>        VectorType;
    typedef CovariantVector<double, 3>    CovariantVectorType;
    typedef vnl_svd< ScalarType >         SolverType;
    typedef vnl_matrix< ScalarType >      InternalMatrixType;
    typedef vnl_vector< ScalarType >      InternalVectorType;

    typedef ProlateSpheroidalTransform<ScalarType> TransformType;

    itkGetObjectMacro (Transform, TransformType);
    itkSetObjectMacro (Transform, TransformType);
  
    itkGetMacro     (UsePiWorkAround, unsigned int);
    itkSetClampMacro(UsePiWorkAround, unsigned int, 0, 1);
    itkBooleanMacro (UsePiWorkAround);

    void GenerateOutputInformation();

    virtual void SetInput( unsigned int, const TMesh* mesh);

  protected:

    ProlateSpheroidalGradientTensorMeshFilter();
    ~ProlateSpheroidalGradientTensorMeshFilter(){};

    void EvaluateAtIndex (unsigned long index);
    void EvaluateUSigmaAnddUl (unsigned long index, InternalMatrixType USigma, InternalMatrixType dUl);
    InternalMatrixType EvaluateSigma (PointType p);
  
    void GenerateData ();
    void ThreadedGenerateData(unsigned long range[2], int threadId );
    void BeforeThreadedGenerateData();
    void AfterThreadedGenerateData();
    int SplitRequestedRegion(int i, int num, unsigned long  range[2]);
    /** Static function used as a "callback" by the MultiThreader.  The threading
     * library will call this routine for each thread, which will delegate the
     * control to ThreadedGenerateData(). */
    static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );

    /** Internal structure used for passing image data into the threading library */
    struct ThreadStruct
    {
      Pointer Filter;
    };
    
    void PrintSelf (std::ostream& os, Indent indent) const
    {
      Superclass::PrintSelf (os, indent);      
    }

  
    InternalVectorType tensor2vec(const TensorType &tensor)
    {
      InternalVectorType vec (TensorType::DegreesOfFreedom);
      for( unsigned int i=0; i<TensorType::DegreesOfFreedom; i++)
	vec[i] = tensor.GetNthComponentAsVector(i);
      return vec;
    }
    TensorType vec2tensor(const InternalVectorType  &vec)
    {
      TensorType tensor (static_cast<ScalarType>(0.0));
      for( unsigned int i=0; i<TensorType::DegreesOfFreedom; i++)
	tensor.SetNthComponentAsVector (i, vec[i]);
      return tensor;
    }
  
    
  private:
    ProlateSpheroidalGradientTensorMeshFilter(const Self&); // purposely not implemented
    void operator=(const Self&); // purposely not implemented
    
    unsigned int               m_UsePiWorkAround;
    typename MeshType::Pointer m_LogInput;
    typename MeshType::Pointer m_LogOutput1;
    typename MeshType::Pointer m_LogOutput2;
    typename MeshType::Pointer m_LogOutput3;

    typename TransformType::Pointer m_Transform;

  };

} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkProlateSpheroidalGradientTensorMeshFilter.txx"
#endif

#endif
