#ifndef _itkDiscreteElasticVectorSmoothingFilter_h_
#define _itkDiscreteElasticVectorSmoothingFilter_h_

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImageToImageFilter.h"

#include "itkAddImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkDiscreteGaussianDerivativeImageFilter.h"
#ifdef ITK_USE_CONCEPT_CHECKING
  #define CONCEPT_CHECKING_IS_USED
  #undef ITK_USE_CONCEPT_CHECKING
#endif
#include "itkDiscreteGaussianImageFilter.h"
#ifdef CONCEPT_CHECKING_IS_USED
  #define ITK_USE_CONCEPT_CHECKING
#endif
#include "itkVectorIndexSelectionCastImageFilter.h"


// So far works only for 3D vectorial images


namespace itk {

template <class TImageType, class TMaskImageType>
class ITK_EXPORT DiscreteElasticVectorSmoothingFilter :
    public ImageToImageFilter<TImageType, TImageType>
{

public:

  /** Usual typedefs */
  typedef DiscreteElasticVectorSmoothingFilter         Self;
  typedef ImageToImageFilter< TImageType, TImageType>  Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;


  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(DiscreteElasticVectorSmoothingFilter, ImageToImageFilter);

  /** Get/Set the standard deviation of the Gaussian smoothing */
  itkGetMacro( Sigma, const double );
  itkSetMacro( Sigma, double );

  /** Get/Set the strength of the elastic filtering */
  itkGetMacro( Kappa, const double );
  itkSetMacro( Kappa, double );

  /** Mask the elastic filtering ON/OFF */
  itkGetMacro( MaskElasticSmoothing, const bool );
  itkSetMacro( MaskElasticSmoothing, bool);
  itkBooleanMacro( MaskElasticSmoothing );

  /** Verbose mode */
  itkGetMacro( Verbose, const bool );
  itkSetMacro( Verbose, bool);
  itkBooleanMacro( Verbose );

  /** Get/Set the image to be used as mask */
  itkSetMacro( MaskImage, typename TMaskImageType::ConstPointer );

  /** Print filter information on the standard output */
  void PrintSelf( std::ostream& os, Indent indent ) const;


protected:

  typedef typename TImageType::PixelType::ValueType  PixelType;

  typedef TMaskImageType MaskImageType;
  typedef typename TMaskImageType::PixelType MaskPixelType;

  typedef  Image< PixelType, 3 > ComponentImageType;

  typedef  VectorIndexSelectionCastImageFilter< TImageType, ComponentImageType >            VectorSplitterType;
  typedef  DiscreteGaussianImageFilter< ComponentImageType, ComponentImageType >            SmootherType;
  typedef  AddImageFilter< ComponentImageType, ComponentImageType, ComponentImageType >     AdderType;
  typedef  DiscreteGaussianDerivativeImageFilter< ComponentImageType, ComponentImageType >  DerivativeFilterType;
  typedef  DivideImageFilter< ComponentImageType, ComponentImageType, ComponentImageType >     DividerType;
  typedef  MultiplyImageFilter< ComponentImageType, ComponentImageType, ComponentImageType >   MultiplierType;
  typedef  MaskImageFilter< ComponentImageType, TMaskImageType, ComponentImageType >        MaskerType;
  typedef  ComposeImageFilter< ComponentImageType >                                 VectorComposerType;

  /// Default constructor
  DiscreteElasticVectorSmoothingFilter();

  /// Generate the data: elastic filtering of the input image
  void GenerateData();


private:

  DiscreteElasticVectorSmoothingFilter(Self&);   // intentionally not implemented
  void operator=(const Self&);           // intentionally not implemented

  /** Standard deviation of the Gaussian smoothing */
  double m_Sigma;

  /** Strength of the elastic filtering */
  double m_Kappa;

  /** Perform the elastic smoothing in a region of interest only */
  bool m_MaskElasticSmoothing;

  /** Display debug messages */
  bool m_Verbose;

  /** Mask image */
  typename MaskImageType::ConstPointer m_MaskImage;
};


}   // end namespace itk




#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiscreteElasticVectorSmoothingFilter.txx"
#endif

#endif
