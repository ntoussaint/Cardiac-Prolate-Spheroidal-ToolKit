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
#ifndef _itk_StackUnstackImageFilter_h_
#define _itk_StackUnstackImageFilter_h_

#include <itkImageToImageFilter.h>

namespace itk
{

  /*! \class StackUnstackImageFilter
    \ingroup ImageAlgorithm
    */

  template<class TInputImage, class TOutputImage>
    class ITK_EXPORT StackUnstackImageFilter :
  public ImageToImageFilter<TInputImage,TOutputImage>
  {
    
  public:
    
    typedef StackUnstackImageFilter                Self;
    typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(StackUnstackImageFilter,ImageToImageFilter);    
    itkStaticConstMacro(OutputImageDimension, unsigned int,
			TOutputImage::ImageDimension);
    itkStaticConstMacro(InputImageDimension, unsigned int,
			TInputImage::ImageDimension);

    
    typedef TInputImage                             InputImageType;
    typedef typename InputImageType::PixelType      InputPixelType;
    typedef typename InputImageType::SpacingType    InputSpacingType;
    typedef typename InputImageType::PointType      InputPointType;
    typedef typename InputImageType::DirectionType  InputDirectionType;
    typedef typename InputImageType::RegionType     InputRegionType;
    typedef typename InputImageType::SizeType       InputSizeType;
    
    typedef TOutputImage                            OutputImageType;
    typedef typename OutputImageType::PixelType     OutputPixelType;
    typedef typename OutputImageType::RegionType    OutputRegionType;
    typedef typename OutputPixelType::ValueType     ScalarType;
    
    
  protected:
    
    StackUnstackImageFilter(){};
    ~StackUnstackImageFilter(){};
    
    void GenerateData(void);    
    /**
       ResampleImageFilter produces an image which is a different size
       than its input.  As such, it needs to provide an implementation
       for GenerateOutputInformation() in order to inform the pipeline
       execution model.  The original documentation of this method is
       below. \sa ProcessObject::GenerateOutputInformaton()
    */
    virtual void GenerateOutputInformation();
    /**
       The current implementation of this class does not supprot
       streaming. As such it produces the output for the largest
       possible region.
    */
    virtual void EnlargeOutputRequestedRegion( DataObject *ptr );

    void PrintSelf(std::ostream& os, Indent indent) const
    {
      Superclass::PrintSelf(os,indent);
    }

  private:
    
    StackUnstackImageFilter(const Self&);
    void operator=(const Self&);

  };



} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStackUnstackImageFilter.txx"
#endif

#endif
