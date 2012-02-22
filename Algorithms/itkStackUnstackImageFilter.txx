/*=========================================================================

Program:   ImagingSciences
Module:    $Id: .h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_StackUnstackImageFilter_txx_
#define _itk_StackUnstackImageFilter_txx_
#include "itkStackUnstackImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <cstdio>

namespace itk
{

  template<class TInputImage, class TOutputImage>
  void
  StackUnstackImageFilter<TInputImage,TOutputImage>
  ::BeforeGenerateData()
  {
  
  }


  template<class TInputImage, class TOutputImage>
  void
  StackUnstackImageFilter<TInputImage,TOutputImage>
  ::GenerateData(void)
  {
  
    typedef ImageRegionIterator<OutputImageType>      IteratorOutputType;
    typedef ImageRegionConstIterator<InputImageType>  IteratorInputType;
    
    unsigned long numPixels = this->GetOutput;
    unsigned long step = numPixels/1000;
    unsigned long progress = 0;
    
    int n = (int)(this->GetNumberOfInputs());
    ScalarType a = 1.0/sqrt(2.0);
    
    IteratorOutputType itOut(this->GetOutput(), outputRegionForThread);

    // create a list of iterators for each input
    std::vector<IteratorInputType> ListOfInputIterators;
    for(int i=0; i<n; i++)
    {
	IteratorInputType it(this->GetInput(i),outputRegionForThread);
	ListOfInputIterators.push_back(it);
    }

    int nonZeroGradientCount = (int)(m_InternalGradientList.size());

    
    if( threadId==0 )
      this->UpdateProgress (0.0);
    
    while(!itOut.IsAtEnd())
    {

      if( this->GetAbortGenerateData() )
        throw itk::ProcessAborted(__FILE__,__LINE__);

      
      OutputPixelType out( static_cast<ScalarType>( 0.0 ) );      
      InputPixelType b0 = ListOfInputIterators[0].Get();
      int nB0 = 1;
 
      // threshold b0
      if(b0>m_BST)
      {
	
        InternalMatrixType B(nonZeroGradientCount, 1);
	int gradCount = 0;

	// any occurence of a null gradient will be considered as a B0
	// and averaged to have an average B0
        for(int i=1; i<n; i++)
        {
	  if (m_GradientList[i-1].GetNorm()>0.001 )
	  {	    
	    ScalarType bi = ListOfInputIterators[i].Get();
	    if( bi<0.001 )
	      bi = 0.001;
	    
	    //B(i-1,0) = log( b0  / bi );
	    B(gradCount++, 0) = bi;
	  }
	  else
	  {  
	    b0 += ListOfInputIterators[i].Get();
	    nB0++;
	  }          
        }

	b0 /= (ScalarType)(nB0);

	for (int i=0; i<nonZeroGradientCount; i++)
	{
	  // it happens that diffusion is greater than the B0: contribution of such gradient is canceled
	  if (b0>B (i,0))
	    B (i,0) = log ( b0 / B (i,0) );
	  else
	    B (i,0) = 0.0;
	}

	
        InternalMatrixType D = m_IG2*m_G*B;

	
        out.SetNthComponent( 0, static_cast<ScalarType>( D(0,0) ));
        out.SetNthComponent( 1, static_cast<ScalarType>( a*D(1,0) ));
        out.SetNthComponent( 2, static_cast<ScalarType>( D(2,0) ));
        out.SetNthComponent( 3, static_cast<ScalarType>( a*D(3,0) ));
        out.SetNthComponent( 4, static_cast<ScalarType>( a*D(4,0) ));
        out.SetNthComponent( 5, static_cast<ScalarType>( D(5,0)) );
	
      }

      if( threadId==0 )
      {
        if( step>0 )
          if( (progress%step)==0 )
            this->UpdateProgress ( double(progress)/double(numPixels) );        
      }
      
      
      itOut.Set(out);
      ++itOut;
      ++progress;
      for(int i=0;i<n;i++)
        ++(ListOfInputIterators[i]);
    }
    if( threadId==0 )
      this->UpdateProgress (1.0);
    
  }


/**
 * Inform pipeline of required output region
 */  
template<class TInputImage, class TOutputImage>
void
StackUnstackImageFilter<TInputImage,TOutputImage>
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // get pointers to the input and output
  ImagePointer outputPtr = const_cast<ImageType*>(this->GetOutput());
  if ( !outputPtr )
  {
    return;
  }

  unsigned int n = this->GetNumberOfInputs();
  if (!n) return;

  unsigned int N1 = this->GetInputImageDimension();
  unsigned int N2 = this->GetOutputImageDimension();
  bool stacking = N1 > N2;
  
  if ( ( !stacking && (n  > 1) ) ||
       (  stacking && (n == 1) ) )
  {
    itkExceptionMacro (<<"Cannot change the image dimension with unconsistent conditions: "
		       <<"N1 = "<<N1"; "
		       <<"N2 = "<<N2"; "
		       <<"n = " <<n);
  }
  
  InputRegionType     inputregion    = this->GetInput (0)->GetLargestPossibleRegion();
  InputSizeType       inputsize      = inputregion.GetSize();
  InputSpacingType    inputspacing   = this->GetInput (0)->GetSpacing();
  InputDirectionType  inputdirection = this->GetInput (0)->GetDirection();
  InputPointType      inputorigin    = this->GetInput (0)->GetOrigin();
  OutputRegionType    outputregion;
  OutputSizeType      outputsize;
  OutputSpacingType   outputspacing;
  OutputDirectionType outputdirection;

   for (unsigned int i=0; i<N2 - 1; i++)
  {
    outputsize[i] = inputsize[i];
    outputspacing[i] = inputspacing[i];
    for (unsigned int j=0; j<N2 - 1; j++)
      outputdirection[i][j] = inputdirection[i][j];
    outputorigin[i] = inputorigin[i];
  }
  
  outputsize[N2 - 1] = n;
  outputregion.SetSize (outputsize);

  if (N1 == N2)
  {
    InputPointType inputorigin2 = this->GetInput (1)->GetOrigin();
    spacing[N2 - 1] = inputorigin2[N1 - 1] - inputorigin[N1 - 1];
    outputdirection = inputdirection;
  }
  
  if (stacking)
  {
    this->SetNumberOfRequiredOutputs (1);

    if (N1 != N2)
      spacing[N2 - 1] = static_cast<ScalarType>(1.0);
    
      
    
    for (unsigned int i=0; i<n; i++)
    {
      if (this->GetInput(i)->GetLargestPossibleRegion()->GetSize() != inputsize)
      {
	itkExceptionMacro (<<"inconsistent size "<<this->GetInput(i)->GetLargestPossibleRegion()->GetSize()
			   <<"on input number "<<i);
      }
    }
  }  
  else
  {
    this->SetNumberOfRequiredOutputs (inputsize[N1 - 1]);

    for(unsigned int i=1; i<inputsize[N1 - 1]; i++)
      this->SetNthOutput (i, OutputImageType::New());
  }
  
  for(unsigned int i=0; i<this->GetNumberOfRequiredOutputs(); i++)
  {
    this->GetOutput (i)->SetLargestPossibleRegion(outputregion);
    this->GetOutput (i)->SetRequestedRegionToLargestPossibleRegion();
  }
  
  // Set spacing and origin
  if (m_UseFixedImage && fixedImage)
  {
    outputPtr->SetSpacing( fixedImage->GetSpacing() );
    outputPtr->SetOrigin( fixedImage->GetOrigin() );
    outputPtr->SetDirection( fixedImage->GetDirection() );
  }
  else
  {
    outputPtr->SetSpacing( movingImage->GetSpacing() );
    outputPtr->SetOrigin( movingImage->GetOrigin() );
    outputPtr->SetDirection( movingImage->GetDirection() );
  }
}

template<class TInputImage, class TOutputImage>
void
StackUnstackImageFilter<TInputImage,TOutputImage>
::EnlargeOutputRequestedRegion(DataObject * ptr )
  {
    // call the superclass's implementation
    Superclass::EnlargeOutputRequestedRegion( ptr );

  }
  
  
} // end of namespace

#endif
