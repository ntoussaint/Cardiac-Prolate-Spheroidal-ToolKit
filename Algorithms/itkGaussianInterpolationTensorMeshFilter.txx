/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkGaussianInterpolationTensorMeshFilter.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_GaussianInterpolationTensorMeshFilter_txx_
#define _itk_GaussianInterpolationTensorMeshFilter_txx_

#include "itkGaussianInterpolationTensorMeshFilter.h"
#include <cstdio>

namespace itk
{

  template<class TMesh>
  void
  GaussianInterpolationTensorMeshFilter<TMesh>
  ::SetInput(unsigned int index, const TMesh* mesh )
  {
    
    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(index, 
                                   const_cast< TMesh *>( mesh ) );
    
  }
  
  template<class TMesh>
  GaussianInterpolationTensorMeshFilter<TMesh>
  ::GaussianInterpolationTensorMeshFilter()
  {
    m_UsePiWorkAround = 0;
    m_Kernel   = GaussianKernelFunctionType::New();
    m_LogInput  = MeshType::New();
    m_LogOutput = MeshType::New();
    m_BandwidthMatrix.SetIdentity();
    m_SqInverseBandwidthMatrix.SetIdentity();
    
    this->SetNumberOfRequiredOutputs(1);
    this->SetNumberOfRequiredInputs(2);

    this->ReleaseDataBeforeUpdateFlagOff();
  };
  
  template<class TMesh>
  void
  GaussianInterpolationTensorMeshFilter<TMesh>
  ::InterpolateAtIndex(unsigned long index)
  {

    MeshType* output = m_LogOutput;

    double sum = 0.0, G = 0.0, deltaX = 0.0;
    PointType pt; pt[0] = pt[1] = pt[2] = 0.0;
    PointType ptn;
    CovariantVectorType dX (static_cast<ScalarType>(0.0));
    TensorType out( static_cast<ScalarType>(0.0) );
    TensorType Mean ( static_cast<ScalarType>(0.0) );
    double epsilon = vcl_numeric_limits<double>::epsilon();
    
    output->GetPoint (index, & pt);
    
    typedef typename MeshType::PointDataContainer  PixelContainer;
    typedef typename MeshType::PointsContainer     PointContainer;
    typename PointContainer::Pointer points  = m_LogInput->GetPoints();
    typename PixelContainer::Pointer tensors = m_LogInput->GetPointData();
    typename PointContainer::Iterator p_it = points->Begin();
    typename PixelContainer::Iterator t_it = tensors->Begin();
    
    while( p_it != points->End())
    {
      ptn = p_it.Value();
      dX = ptn - pt;
      
      if (m_UsePiWorkAround)
      {
	while (dX[2] >= vnl_math::pi)
	  dX[2] -= 2 * vnl_math::pi;
	while (dX[2] < - vnl_math::pi)
	  dX[2] += 2 * vnl_math::pi;
	while (dX[1] >= vnl_math::pi)
	  dX[1] -= 2 * vnl_math::pi;
	while (dX[1] < - vnl_math::pi)
	  dX[1] += 2 * vnl_math::pi;
      }

      deltaX = std::sqrt ( dX * ( m_SqInverseBandwidthMatrix * dX ) );
      
      G = this->GetKernel()->Evaluate (deltaX);
      //G = 1;
      
      sum += G;
      Mean += t_it.Value() * G ;

      ++p_it;
      ++t_it;
    }
    
    // We have to be careful here that the mean tensor is not too close
    // to the null tensor, and that sum is greater than zero.
    
    if (sum > epsilon)
    {
      Mean /= static_cast<ScalarType>(sum);
      out = Mean;
      output->SetPointData (index, out);
    }
    else
    {
      std::ostringstream os;
      os<<"fj was too small : "<<"sum of weights is null";
      itkWarningMacro (<<os.str().c_str());      
    }

    
  }
  
  
  template<class TMesh>
  void
  GaussianInterpolationTensorMeshFilter<TMesh>
  ::BeforeThreadedGenerateData()
  {
    std::cout<<"tensors LOG"<<std::endl;
    
    typedef typename MeshType::PointDataContainer  PixelContainer;
    typename PixelContainer::ConstPointer pixels = this->GetInput(1)->GetPointData();
    typename PixelContainer::Pointer logpixels = m_LogInput->GetPointData();
    typename PixelContainer::ConstIterator it = pixels->Begin();
    typename PixelContainer::Iterator log_it = logpixels->Begin();
    
    while( it != pixels->End() )
    {
      if (!it.Value().IsFinite() || it.Value().HasNans())
	std::cout<<"T is given not finite at "<<it.Value()<<std::endl;
      else
	log_it.Value() = it.Value().Log();
      
      ++it;
      ++log_it;
    }
  }
  
  template<class TMesh>
  void GaussianInterpolationTensorMeshFilter<TMesh>
  ::AfterThreadedGenerateData()
  {

    std::cout<<"tensors EXP"<<std::endl;

    typedef typename MeshType::PointDataContainer  PixelContainer;
    typename PixelContainer::Pointer pixels    = this->GetOutput()->GetPointData();
    typename PixelContainer::Pointer logpixels =       m_LogOutput->GetPointData();
    typename PixelContainer::Iterator it    =    pixels->Begin();
    typename PixelContainer::Iterator logit = logpixels->Begin();

    while( logit != logpixels->End() )
    {
      if (!logit.Value().IsFinite() || logit.Value().HasNans())
	std::cout<<"T is given not finite at "<<logit.Value()<<std::endl;
      else
	it.Value() = logit.Value().Exp();
      
      ++it; ++logit;
      
    }
    
  }

  
  template<class TMesh>
  void 
  GaussianInterpolationTensorMeshFilter<TMesh>
  ::GenerateData()
  {
    std::cout<<"Gaussian interpolation starting"<<std::endl;
    std::cout<<"using bandwidth matrix \n"<<m_BandwidthMatrix<<std::endl;
    
    // Call a method that can be overridden by a subclass to perform
    // some calculations prior to splitting the main computations into
    // separate threads
    this->BeforeThreadedGenerateData();
    
    // Set up the multithreaded processing
    ThreadStruct str;
    str.Filter = this;
    
    this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
    this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str);

    std::cout<<"interpolation..."<<std::endl;
    std::cout<<"|";
    
    // multithread the execution
    this->GetMultiThreader()->SingleMethodExecute();

    std::cout<<"|"<<std::endl<<"done."<<std::endl;

    // Call a method that can be overridden by a subclass to perform
    // some calculations after all the threads have completed
    this->AfterThreadedGenerateData();

  }
  
  template<class TMesh>
  void
  GaussianInterpolationTensorMeshFilter<TMesh>
  ::GenerateOutputInformation()
  {
    typedef typename MeshType::PointsContainer    PointContainer;
    typedef typename MeshType::PointDataContainer PixelContainer;
    
    // call the superclass's implementation of this method
    Superclass::GenerateOutputInformation();
    
    typename MeshType::Pointer      output  = this->GetOutput();
    typename MeshType::ConstPointer input1  = this->GetInput(0);
    typename MeshType::ConstPointer input2  = this->GetInput(1);

    if( !input1)
    {
      itkExceptionMacro(<<"Missing Input Mesh (sampling points)");
    }
    
    if( !input2)
    {
      itkExceptionMacro(<<"Missing Input Mesh (data)");
    }
    
    if( !output )
    {
      itkExceptionMacro(<<"Missing Output Mesh");
    }

    this->CopyInputMeshToOutputMeshPoints();
    
    typename PixelContainer::Pointer outputdata   = PixelContainer::New();    
    outputdata->Reserve (input1->GetNumberOfPoints());
    output->SetPointData (outputdata);

    typename PointContainer::ConstPointer datapoints = input2->GetPoints();
    typename PointContainer::ConstPointer samplingpoints = input1->GetPoints();
    typename PixelContainer::ConstPointer datapixels = input2->GetPointData();
    typename PointContainer::Pointer logdatapoints1 = PointContainer::New();
    typename PointContainer::Pointer logdatapoints2 = PointContainer::New();
    typename PixelContainer::Pointer logdatapixels1 = PixelContainer::New();
    typename PixelContainer::Pointer logdatapixels2 = PixelContainer::New();
    logdatapoints1->Reserve (datapoints->Size());
    logdatapoints2->Reserve (samplingpoints->Size());
    logdatapixels1->Reserve (datapixels->Size());
    logdatapixels2->Reserve (samplingpoints->Size());
    
    m_LogInput->SetPointData  (logdatapixels1);
    m_LogOutput->SetPointData (logdatapixels2);
    m_LogInput->SetPoints  (logdatapoints1);
    m_LogOutput->SetPoints (logdatapoints2);

    typename PointContainer::ConstIterator  in_it  = datapoints->Begin();
    typename PointContainer::Iterator      out_it  = logdatapoints1->Begin();
    while(in_it != datapoints->End())
    {
      out_it.Value() = in_it.Value();
      ++in_it; ++out_it;
    }
    in_it  = samplingpoints->Begin();
    out_it = logdatapoints2->Begin();
    while(in_it != samplingpoints->End())
    {
      out_it.Value() = in_it.Value();
      ++in_it; ++out_it;
    }
    
  }
  
  template<class TMesh>
  void
  GaussianInterpolationTensorMeshFilter<TMesh>
  ::ThreadedGenerateData(unsigned long range[2], int id)
  {
    unsigned long start = range[0];
    unsigned long end   = start + range[1];

    double percent = 0;
    int done = 0;
    double threads = this->GetNumberOfThreads();
    
    for (unsigned long i = start; i < end; i++)
    {

      percent += 100 * 1.0 / (double)( threads * (end - start) );
      if (std::floor (percent) >= done)
      {
	done++;
	std::cout<<"="<<std::flush;
      }
      try
      {
	this->InterpolateAtIndex (i);
      }
      catch (itk::ExceptionObject & e)
      {
	std::cerr << e << std::endl;
	continue;
      }
    }
  }

  template<class TMesh>
  int
  GaussianInterpolationTensorMeshFilter<TMesh>
  ::SplitRequestedRegion(int i, int num, unsigned long range[2])
  {
    MeshType* output  = this->GetOutput();

    if( !output )
    {
      itkExceptionMacro(<<"Missing Output Mesh");
    }
    
    // determine the actual number of pieces that will be generated
    unsigned long totalnumber;
    totalnumber = output->GetNumberOfPoints();
    
    int valuesPerThread = (int)::vcl_ceil(totalnumber/(double)num);
    int maxThreadIdUsed = (int)::vcl_ceil(totalnumber/(double)valuesPerThread) - 1;

    range[0] = 0;
    range[1] = totalnumber;
    
    if (i < maxThreadIdUsed)
    {
      range[0] += i*valuesPerThread;
      range[1] = valuesPerThread;
    }
    if (i == maxThreadIdUsed)
    {
      range[0] += i*valuesPerThread;
      range[1] = range[1] - i*valuesPerThread;
    }
  
    return maxThreadIdUsed + 1;
  }  


  // Callback routine used by the threading library. This routine just calls
  // the ThreadedGenerateData method after setting the correct region for this
  // thread. 
  
  template<class TMesh>
  ITK_THREAD_RETURN_TYPE
  GaussianInterpolationTensorMeshFilter<TMesh>
  ::ThreaderCallback( void *arg )
  {
    ThreadStruct *str;
    int total, threadId, threadCount;
    
    threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
    threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;
    
    str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);
    
    // execute the actual method with appropriate output region
    // first find out how many pieces extent can be split into.
    unsigned long range[2];
    
    total = str->Filter->SplitRequestedRegion(threadId, threadCount, range);
    
    if (threadId < total)
    {
      str->Filter->ThreadedGenerateData(range, threadId);
    }
    // else
    //   {
    //   otherwise don't use this thread. Sometimes the threads dont
    //   break up very well and it is just as efficient to leave a 
    //   few threads idle.
    //   }
    
    return ITK_THREAD_RETURN_VALUE;
  }
  
  
} // end of namespace itk

#endif
