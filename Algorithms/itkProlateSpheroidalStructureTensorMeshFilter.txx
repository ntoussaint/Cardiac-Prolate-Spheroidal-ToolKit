/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkProlateSpheroidalStructureTensorMeshFilter.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_ProlateSpheroidalStructureTensorMeshFilter_txx_
#define _itk_ProlateSpheroidalStructureTensorMeshFilter_txx_

#include "itkProlateSpheroidalStructureTensorMeshFilter.h"
#include <cstdio>

namespace itk
{

  template<class TMesh>
  void
  ProlateSpheroidalStructureTensorMeshFilter<TMesh>
  ::SetInput(unsigned int index, const TMesh* mesh )
  {
    
    // Process object is not const-correct so the const_cast is required here
    this->ProcessObject::SetNthInput(index, 
                                   const_cast< TMesh *>( mesh ) );
    
  }
  
  template<class TMesh>
  ProlateSpheroidalStructureTensorMeshFilter<TMesh>
  ::ProlateSpheroidalStructureTensorMeshFilter()
  {
    m_UsePiWorkAround = 1;
    m_Transform = 0;

    m_GradientFilter = GradientFilterType::New();
    
    this->SetNumberOfRequiredOutputs(1);
    this->SetNumberOfRequiredInputs(2);

    this->ReleaseDataBeforeUpdateFlagOff();
  };
  
  template<class TMesh>
  void
  ProlateSpheroidalStructureTensorMeshFilter<TMesh>
  ::EvaluateAtIndex(unsigned long index)
  {
    typename MeshType::Pointer      output  = this->GetOutput();
    
    TensorType g1 (static_cast<ScalarType>(0.0));
    TensorType g2 (static_cast<ScalarType>(0.0));
    TensorType g3 (static_cast<ScalarType>(0.0));
    
    m_Gradient1->GetPointData (index, &g1);
    m_Gradient2->GetPointData (index, &g2);
    m_Gradient3->GetPointData (index, &g3);

    InternalMatrixType gradl (3, TensorType::DegreesOfFreedom);
    gradl.set_row (0, this->tensor2vec (g1));
    gradl.set_row (1, this->tensor2vec (g2));
    gradl.set_row (2, this->tensor2vec (g3));

    InternalMatrixType T = gradl * gradl.transpose();

    TensorType out;
    out.SetVnlMatrix (T);
    
    if (out.GetFA() > 0.95)
    {
      itkWarningMacro (<<"outlier tensor : \n"<<out);
      out = static_cast<TensorType>(0.0);
    }
    
    output->SetPointData (index, out);
  }

  template<class TMesh>
  void
  ProlateSpheroidalStructureTensorMeshFilter<TMesh>
  ::BeforeThreadedGenerateData()
  {

    m_GradientFilter->SetInput (0, this->GetInput (0));
    m_GradientFilter->SetInput (1, this->GetInput (1));
    m_GradientFilter->SetTransform (this->GetTransform());
    m_GradientFilter->SetUsePiWorkAround (this->GetUsePiWorkAround());

    m_GradientFilter->Update();

    m_Gradient1 = m_GradientFilter->GetOutput (0);
    m_Gradient2 = m_GradientFilter->GetOutput (1);
    m_Gradient3 = m_GradientFilter->GetOutput (2);

    m_Gradient1->DisconnectPipeline();
    m_Gradient2->DisconnectPipeline();
    m_Gradient3->DisconnectPipeline();
    
    std::cout<<"structure: tensors LOG"<<std::endl;
    
    typedef typename MeshType::PointDataContainer  PixelContainer;
    typename PixelContainer::Pointer pixels1    = m_Gradient1->GetPointData();
    typename PixelContainer::Pointer pixels2    = m_Gradient2->GetPointData();
    typename PixelContainer::Pointer pixels3    = m_Gradient3->GetPointData();
    typename PixelContainer::Iterator it1       = pixels1->Begin();
    typename PixelContainer::Iterator it2       = pixels2->Begin();
    typename PixelContainer::Iterator it3       = pixels3->Begin();

    while( it1 != pixels1->End() )
    {
      if (!it1.Value().IsFinite() || it1.Value().HasNans() || !it1.Value().IsPositive())
	std::cout<<"T1 is given not finite at "<<it1.Value()<<std::endl;
      else
	it1.Value() = it1.Value().Log();
      if (!it2.Value().IsFinite() || it2.Value().HasNans() || !it2.Value().IsPositive())
	std::cout<<"T2 is given not finite at "<<it2.Value()<<std::endl;
      else
	it2.Value() = it2.Value().Log();
      if (!it3.Value().IsFinite() || it3.Value().HasNans() || !it3.Value().IsPositive())
	std::cout<<"T3 is given not finite at "<<it3.Value()<<std::endl;
      else
	it3.Value() = it3.Value().Log();
      ++it1; ++it2; ++it3;
    }
    
  }
  
  template<class TMesh>
  void ProlateSpheroidalStructureTensorMeshFilter<TMesh>
  ::AfterThreadedGenerateData()
  {

    std::cout<<"structure: tensors EXP"<<std::endl;

    typedef typename MeshType::PointDataContainer  PixelContainer;
    typename PixelContainer::Pointer pixels    = this->GetOutput()->GetPointData();
    typename PixelContainer::Iterator it       = pixels->Begin();

    while( it != pixels->End() )
    {
      if (!it.Value().IsFinite() || it.Value().HasNans())
	std::cout<<"T is given not finite at "<<it.Value()<<std::endl;
      else
	it.Value() = it.Value().Exp();
      ++it;
    }
    
  }

  
  template<class TMesh>
  void 
  ProlateSpheroidalStructureTensorMeshFilter<TMesh>
  ::GenerateData()
  {
    std::cout<<"Structure Evaluation starting"<<std::endl;
    
    // Call a method that can be overridden by a subclass to perform
    // some calculations prior to splitting the main computations into
    // separate threads
    this->BeforeThreadedGenerateData();
    
    // Set up the multithreaded processing
    ThreadStruct str;
    str.Filter = this;
    
    this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
    this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str);

    std::cout<<"structure-evaluation..."<<std::endl;
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
  ProlateSpheroidalStructureTensorMeshFilter<TMesh>
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
    
    if( !output)
    {
      itkExceptionMacro(<<"Missing Output Meshes");
    }

    this->CopyInputMeshToOutputMeshPoints();
    
    typename PixelContainer::Pointer outputdata   = PixelContainer::New();    
    outputdata->Reserve (input1->GetNumberOfPoints());
    output->SetPointData (outputdata);
    
  }
  
  template<class TMesh>
  void
  ProlateSpheroidalStructureTensorMeshFilter<TMesh>
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
	this->EvaluateAtIndex (i);
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
  ProlateSpheroidalStructureTensorMeshFilter<TMesh>
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
  ProlateSpheroidalStructureTensorMeshFilter<TMesh>
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
