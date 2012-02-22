#ifndef _itk_EigenVectorToTensorImageFilter_h_
#define _itk_EigenVectorToTensorImageFilter_h_

#include <itkImageToImageFilter.h>

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT EigenVectorToTensorImageFilter :
public ImageToImageFilter<TInputImage, TOutputImage>
{

public:

    typedef EigenVectorToTensorImageFilter                      Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self>                            Pointer;
    typedef SmartPointer<const Self>                      ConstPointer;

    itkNewMacro(  Self );
    itkTypeMacro( EigenVectorToTensorImageFilter, ImageToImageFilter );

    typedef TInputImage                                   InputImageType;
    typedef typename InputImageType::PixelType            InputPixelType;
    typedef TOutputImage                                  OutputImageType;
    typedef typename OutputImageType::PixelType           OutputPixelType;
    typedef typename OutputImageType::RegionType          OutputImageRegionType;


protected:

    EigenVectorToTensorImageFilter(void){};
    virtual ~EigenVectorToTensorImageFilter(void){};

    void BeforeThreadedGenerateData(void);
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int threadId);
    void PrintSelf(std::ostream& os, Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }


private:


};


} // end of namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEigenVectorToTensorImageFilter.txx"
#endif

#endif
