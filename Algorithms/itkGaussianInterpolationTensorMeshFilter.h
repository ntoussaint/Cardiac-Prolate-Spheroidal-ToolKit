/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkGaussianInterpolationTensorMeshFilter.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_GaussianInterpolationTensorMeshFilter_h_
#define _itk_GaussianInterpolationTensorMeshFilter_h_

#include <itkMeshToMeshFilter.h>
#include <itkMesh.h>
#include <itkTensor.h>
#include <itkPoint.h>
#include <itkGaussianKernelFunction.h>

#include <vector>

namespace itk
{

  /**
     \class GaussianInterpolationTensorMeshFilter
     \brief This class is a kernel based interpolation for tensors embedded in a mesh grid

     GaussianInterpolationTensorMeshFilter is a mesh to mesh filter that takes two inputs.
     Use SetInput(0, ...) to give the data to be used for interpolation. Use SetInput(1, ...)
     to give the mesh points where the interpolation is performed. The output mesh has the same
     points and structure as the second input (sampling points). At each point position of the
     output mesh, an anisotropic kernel based mean of all tensors given in the first input is
     estimated. 

     The kernel estimate can be written as followed:
     \f[
     \hat{m}_{\sigma}(D_{\xi}) = \exp \left( \frac{ \sum_{i=1}^N K_{\sigma} (\xi - \xi_i) \log(D_{\xi_i})}{\sum_{i=1}^N
     K_{\sigma}(\xi -\xi_i)} \right)
     \f]

     We interpolate tensors in the tensor log space.

     Use SetAlpha() to control the anisotropy of the kernel in the three degrees of freedom. 

     In case of Non-Cartesian coordinates spatial positions (spherical, cylindrical, or prolate
     spheroidal), you should use UsePiWorkAroundOn(). This method will trigger a check during
     the spatial distance estimation, and module around -pi and pi for the second and third component,
     as the first one is supposed to be a radial distance.
     
     \sa GaussianKernelFunction MeshToMeshFilter
     Author: Nicolas Toussaint. Copyright KCL 2010.
  */
  
  template <class TMesh = Mesh < Tensor<double,3>, 3> >
  class ITK_EXPORT GaussianInterpolationTensorMeshFilter :
  public MeshToMeshFilter< TMesh, TMesh>
  {
    
  public:

    typedef GaussianInterpolationTensorMeshFilter Self;
    typedef TMesh MeshType;
    typedef MeshToMeshFilter< MeshType, MeshType > Superclass;

    typedef typename MeshType::Pointer    MeshPointer;
    typedef typename MeshType::PixelType  PixelType;

    typedef SmartPointer<Self>         Pointer;
    typedef SmartPointer<const Self>   ConstPointer;
    
    itkNewMacro(Self);
    itkTypeMacro(GaussianInterpolationTensorMeshFilter, Superclass);
 
    /** Image typedefs */
    typedef PixelType TensorType;
    typedef typename PixelType::ValueType ScalarType;
    typedef typename MeshType::PointType   PointType;
    typedef Vector <ScalarType, 3>   VectorType;
    typedef KernelFunction KernelFunctionType;
    typedef typename KernelFunctionType::Pointer KernelFunctionPointerType;
    typedef GaussianKernelFunction GaussianKernelFunctionType;
    typedef CovariantVector<double, 3> CovariantVectorType;
    typedef Matrix<double, 3, 3> BandwidthMatrixType;
  
    void SetAlpha (double alpha[3])
    {
      
      // the 3 parameters correspond to the eigen values
      // of the diagonal bandwidth matrix H.
      m_BandwidthMatrix.SetIdentity();
      
      for (unsigned int i=0; i<3; i++)
    	m_BandwidthMatrix[i][i] = alpha[i];
      // the inverse bandwidth matrix
      // is estimated.
      m_SqInverseBandwidthMatrix = m_BandwidthMatrix.GetInverse ();
      m_SqInverseBandwidthMatrix *= m_SqInverseBandwidthMatrix;
      this->Modified();
    }
  
    itkGetMacro (BandwidthMatrix, BandwidthMatrixType);
    void SetBandwidthMatrix (BandwidthMatrixType matrix)
    {
      m_BandwidthMatrix = matrix;
      // the inverse bandwidth matrix
      // is estimated.
      m_SqInverseBandwidthMatrix = m_BandwidthMatrix.GetInverse ();
      m_SqInverseBandwidthMatrix *= m_SqInverseBandwidthMatrix;
      this->Modified();
    }
    
    itkGetMacro     (UsePiWorkAround, unsigned int);
    itkSetClampMacro(UsePiWorkAround, unsigned int, 0, 1);
    itkBooleanMacro (UsePiWorkAround);

    itkGetConstObjectMacro (Kernel, KernelFunctionType);
    itkSetObjectMacro (Kernel, KernelFunctionType);

    void GenerateOutputInformation();

    virtual void SetInput( unsigned int, const TMesh* mesh);

  protected:

    GaussianInterpolationTensorMeshFilter();
    ~GaussianInterpolationTensorMeshFilter(){};

    void InterpolateAtIndex (unsigned long index);
    
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
    
    
  private:
    GaussianInterpolationTensorMeshFilter(const Self&); // purposely not implemented
    void operator=(const Self&); // purposely not implemented
    
    unsigned int               m_UsePiWorkAround;
    BandwidthMatrixType        m_BandwidthMatrix;
    BandwidthMatrixType        m_SqInverseBandwidthMatrix;
    KernelFunctionPointerType  m_Kernel;
    typename MeshType::Pointer m_LogInput;
    typename MeshType::Pointer m_LogOutput;
    

  };

} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaussianInterpolationTensorMeshFilter.txx"
#endif

#endif
