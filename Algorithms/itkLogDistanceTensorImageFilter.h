/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkLogDistanceTensorImageFilter.h 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_LogDistanceTensorImageFilter_h_
#define _itk_LogDistanceTensorImageFilter_h_


#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkTensor.h"

namespace itk
{
  /*! \class LogDistanceTensorImageFilter
    \ingroup 
  */

  template<class TInputImage, class TOutputImage>
    class ITK_EXPORT LogDistanceTensorImageFilter :
  public ImageToImageFilter< TInputImage, TOutputImage >
  {

  public:
    /** Standard class typedefs */
    typedef LogDistanceTensorImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    itkNewMacro (Self);
    itkTypeMacro (LogDistanceTensorImageFilter, ImageToImageFilter);

    itkStaticConstMacro (ImageDimension, unsigned int, TOutputImage::ImageDimension);

    /** Image typedefs */
    typedef TInputImage                        InputImageType;
    typedef typename TInputImage::Pointer      InputImagePointer;
    typedef typename TInputImage::ConstPointer InputImageConstPointer;
    typedef TOutputImage                        OutputImageType;
    typedef typename TOutputImage::Pointer      OutputImagePointer;
    typedef typename TOutputImage::ConstPointer OutputImageConstPointer;
    
    typedef typename TOutputImage::PixelType        OutputPixelType;
    typedef typename TOutputImage::PointType        PointType;
    typedef typename PointType::VectorType          VectorType;
    typedef typename TInputImage::PixelType         InputPixelType;
    typedef typename InputPixelType::RealValueType  RealType;
    typedef Image <RealType, ImageDimension> DistanceImageType;    
    typedef typename DistanceImageType::Pointer DistanceImagePointerType;    
    typedef typename DistanceImageType::PixelType DistancePixelType;    
    
    /** Image typedefs support */
    
    /** Superclass typedefs */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;


    
    /** Set the first input. */
    void SetInput1( const InputImageType * image )
    { this->SetInput( image ); }

    /** Set the second input. */
    void SetInput2( const InputImageType * image );

    /** Get the first input. */
    const InputImageType * GetInput1(void)
    { return this->GetInput(); }
  
    /** Get the second input. */
    const InputImageType * GetInput2(void);

    /** Get the output. */
    const DistanceImageType * GetOutput(void)
    { return this->GetDistanceImage(); }
    
    virtual void GenerateInputRequestedRegion() throw (InvalidRequestedRegionError);

    itkGetConstObjectMacro( DistanceImage, DistanceImageType );
    
  protected:

    LogDistanceTensorImageFilter();
    virtual ~LogDistanceTensorImageFilter(){}
    
    void ThreadedGenerateData (const OutputImageRegionType&, int);
    void PrintSelf (std::ostream&, Indent) const;
    
  private:
    LogDistanceTensorImageFilter (const Self&);
    void operator=(const Self&);

    DistanceImagePointerType m_DistanceImage;
    
    
  };
  


} // end of namespace


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLogDistanceTensorImageFilter.txx"
#endif

#endif
