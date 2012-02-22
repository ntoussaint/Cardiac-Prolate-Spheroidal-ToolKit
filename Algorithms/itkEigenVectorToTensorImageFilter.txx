#ifndef _itk_EigenVectorToTensorImageFilter_txx_
#define _itk_EigenVectorToTensorImageFilter_txx_

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include "itkEigenVectorToTensorImageFilter.h"


namespace itk
{


template <class TInputImage, class TOutputImage>
void
EigenVectorToTensorImageFilter<TInputImage,TOutputImage>::
BeforeThreadedGenerateData()
{
    // Preliminary test
    if( InputPixelType::Dimension != 3 )
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error: Tensor and vector dimensions do not match.");
    if( OutputPixelType::Dimension != 3 )
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error: Tensor and vector dimensions do not match.");
  
    // Images
    const TInputImage * input  = this->GetInput();
    TOutputImage      * output = this->GetOutput();

    // Copy image geometry from input to output image
    output->SetRegions(   input->GetLargestPossibleRegion() );
    output->SetOrigin(    input->GetOrigin() );
    output->SetSpacing(   input->GetSpacing() );
    output->SetDirection( input->GetDirection() );
}



template <class TInputImage, class TOutputImage>
void
EigenVectorToTensorImageFilter<TInputImage,TOutputImage>::
ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int threadId)
{

    // Define iterators
    typedef ImageRegionConstIterator<InputImageType>  IteratorInputType;
    typedef ImageRegionIterator<OutputImageType>      IteratorOutputType;
    IteratorInputType  itIn (this->GetInput(),  outputRegionForThread);
    IteratorOutputType itOut(this->GetOutput(), outputRegionForThread);
    typedef typename OutputPixelType::ValueType ValueType;
    typedef CrossHelper<InputPixelType> CrossType;
    CrossType crosser;
    InputPixelType  V1, V2, V3, Ox (0.0), Oy (0.0);
    Ox[0] = 1; Oy[1] = 1.0;
    
    V1[0] = 0; V1[1] = 1; V1[2] = 0;
    
    ValueType e1, e2, e3;
    OutputPixelType T;
    const ValueType epsilon = static_cast<ValueType>(0.01);
    
    // Iterate on image voxels
    while( !itOut.IsAtEnd() )
    {
      T = static_cast<OutputPixelType>(0.0);
      
      if (itIn.Get().GetNorm() > epsilon)
      {
	
        V1 = itIn.Get();
        e1 = V1.GetNorm(); V1.Normalize();
	if ( ( V1 * Ox ) >= (1.0 - epsilon) )
	  V2 = crosser (V1, Oy);
	else
	  V2 = crosser (V1, Ox);
	
	V3 = crosser (V1, V2);
	e2 = 0.6 * e1;
	e3 = 0.6 * e1;

	// std::cout<<"V1: "<<V1<<std::endl;
	// std::cout<<"V1 * Ox : "<<(double)(V1 * Ox)<<std::endl;
	// std::cout<<"V2: "<<V2<<std::endl;
	// std::cout<<"V3: "<<V3<<std::endl;
	
	T = OutputPixelType (e1 * V1) + OutputPixelType (e2 * V2) + OutputPixelType (e3 * V3);

	// std::cout<<std::endl<<"T: "<<std::endl<<T<<std::endl;
	// getchar();
      }
      
      itOut.Set(T);
      
      ++itOut;
      ++itIn;
    }
    
}


} // end of namespace


#endif
