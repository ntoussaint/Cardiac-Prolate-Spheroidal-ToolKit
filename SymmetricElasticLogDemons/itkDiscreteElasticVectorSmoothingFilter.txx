#ifndef _itkDiscreteElasticVectorSmoothingFilter_txx
#define _itkDiscreteElasticVectorSmoothingFilter_txx

#include "itkDiscreteElasticVectorSmoothingFilter.h"

#include "itkAddImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkCompose3DVectorImageFilter.h"
#include "itkDivideByConstantImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDiscreteGaussianDerivativeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkDerivativeImageFilter.h"



namespace itk
{

// -----------------------------------------------------------------------------
// Constructor

template < class TImageType, class TMaskImageType >
DiscreteElasticVectorSmoothingFilter< TImageType, TMaskImageType >::
DiscreteElasticVectorSmoothingFilter()
{
  m_Sigma = 1.0;
  m_Kappa = 0.5;
  m_Verbose = false;
  m_MaskElasticSmoothing = false;
  m_MaskImage = NULL;
}



// -----------------------------------------------------------------------------
// GenerateData

template < class TImageType, class TMaskImageType >
void
DiscreteElasticVectorSmoothingFilter< TImageType, TMaskImageType >::
GenerateData()
{

  // Extract the components of the deformation field
  typename ComponentImageType::Pointer F[3];
  for ( unsigned int i = 0; i < 3; ++i )
  {
    typename VectorSplitterType::Pointer extractor = VectorSplitterType::New();
    extractor->SetInput( this->GetInput() );
    extractor->SetIndex( i );
    extractor->Update();
    F[i] = extractor->GetOutput();
  }


  if ( m_Verbose )
    std::cerr << "+ Gaussian smoothing, sigma=" << m_Sigma << "\n";

  // Initial Gaussian smoothing
  typename ComponentImageType::Pointer S[3];
  for ( unsigned int i = 0; i < 3; ++i )
  {
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetInput( F[i] );
    smoother->SetVariance( m_Sigma*m_Sigma );
//    smoother->SetMaximumError( 0.1 );
//    smoother->SetMaximumKernelWidth( 100 );
    smoother->SetUseImageSpacingOn();
    smoother->Update();
    S[i] = smoother->GetOutput();
  }

  // Elastic smoothing
  if ( m_Kappa != 0.0 )
  {
    if ( m_Verbose )
      std::cerr << "+ Elastic smoothing, kappa=" << m_Kappa << "\n";

    typename ComponentImageType::SpacingType spacing = this->GetInput()->GetSpacing();

    // Caches used during the computation of the elastic smoothing
    typename ComponentImageType::Pointer H[3];
    for ( unsigned int i = 0; i < 3; ++i )
    {
      H[i]= ComponentImageType::New();
      H[i]->CopyInformation( F[i] );
      H[i]->SetRegions( F[i]->GetRequestedRegion() );
      H[i]->Allocate();
      H[i]->FillBuffer( 0.0 );
    }

    unsigned int orders[3];

    for ( unsigned int i = 0; i < 3; ++i )   // Rows of the Hessian
    {
      for ( unsigned int j = 0; j < 3; ++j ) // Columns of the Hessian
      {
        // Get the derivation order of each direction according to the
        // matrix position
        orders[0] = 0;
        orders[1] = 0;
        orders[2] = 0;

        orders[i] = orders[i] + 1;
        orders[j] = orders[j] + 1;

        // Differentiate
        typename ComponentImageType::Pointer diff = S[j];
        for ( unsigned int o = 0; o < 3; ++o )
        {
          typedef DerivativeImageFilter< ComponentImageType, ComponentImageType > DerivatorFilterType;
          typename DerivatorFilterType::Pointer derivator = DerivatorFilterType::New();
          derivator->SetInput( diff );
          derivator->SetOrder( orders[o] );
          derivator->SetDirection( o );
          derivator->SetUseImageSpacingOn();
          derivator->Update();

          diff = derivator->GetOutput();
          diff->DisconnectPipeline();
        }

        // Adding
        typename AdderType::Pointer adder = AdderType::New();
        adder->SetInput1( H[i] );
        adder->SetInput2( diff );
        adder->InPlaceOn();
        adder->Update();

        H[i] = adder->GetOutput();
        H[i]->DisconnectPipeline();

      } // End loop on the columns j
    } // End loop on the lines i

    // Mask the elastic component
    if ( m_MaskElasticSmoothing && m_MaskImage )
    {
      if ( m_Verbose )
        std::cerr << "+ Masking\n";
      for ( unsigned int i = 0; i < 3; ++i )
      {
        typename MaskerType::Pointer masker = MaskerType::New();
        masker->SetInput1( H[i] );
        masker->SetInput2( m_MaskImage );
        masker->Update();
        H[i] = masker->GetOutput();
      }
    }

    // Add the resulting elastic component region to the Gaussian filter
    for ( unsigned int i = 0; i < 3; ++i )
    {
      typename MultiplierType::Pointer multiplier = MultiplierType::New();
      multiplier->SetInput( H[i] );
      multiplier->SetConstant( m_Sigma * m_Sigma * m_Kappa / (1.0 + m_Kappa) );

      typename AdderType::Pointer adder = AdderType::New();
      adder->InPlaceOn();
      adder->SetInput1( S[i] );
      adder->SetInput2( multiplier->GetOutput() );
      adder->Update();
      S[i] = adder->GetOutput();
      S[i]->DisconnectPipeline();
    }

  } // End elastic smoothing

  // Rebuild the vector field
  typename VectorComposerType::Pointer ComposeVectorFilter = VectorComposerType::New();
  ComposeVectorFilter->SetInput1( S[0] );
  ComposeVectorFilter->SetInput2( S[1] );
  ComposeVectorFilter->SetInput3( S[2] );
  ComposeVectorFilter->GraftOutput( this->GetOutput() );
  ComposeVectorFilter->Update();

  this->GraftOutput( ComposeVectorFilter->GetOutput() );
  this->GetOutput()->Modified();
}



// -----------------------------------------------------------------------------
// PrintSelf method

template <class TImageType, class TMaskImageType >
void
DiscreteElasticVectorSmoothingFilter< TImageType, TMaskImageType >::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Sigma: "
      << m_Sigma << std::endl;

  os << indent << "Kappa: "
      << m_Kappa << std::endl;

  os << indent << "Mask elastic smoothing: "
      << m_MaskElasticSmoothing << std::endl;

  os << indent << "Verbose: "
      << m_Verbose << std::endl;
}



} // end namespace itk

#endif

