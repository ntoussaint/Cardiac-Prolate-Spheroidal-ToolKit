#ifndef _itkElasticVectorSmoothingFilter_txx
#define _itkElasticVectorSmoothingFilter_txx

#include "itkElasticVectorSmoothingFilter.h"

#include "itkAddImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkCompose3DVectorImageFilter.h"
#include "itkDivideByConstantImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"



namespace itk
{

// -----------------------------------------------------------------------------
// Constructor

template < class TImageType, class TMaskImageType >
ElasticVectorSmoothingFilter< TImageType, TMaskImageType >::
ElasticVectorSmoothingFilter()
{
  m_Sigma = 1.0;
  m_Kappa = 0.5;
  m_Verbose = false;
  m_MaskElasticSmoothing = false;
  m_MaskImage = NULL;
  m_NormalizeAcrossScale = false;
}



// -----------------------------------------------------------------------------
// GenerateData

template < class TImageType, class TMaskImageType >
void
ElasticVectorSmoothingFilter< TImageType, TMaskImageType >::
GenerateData()
{

  std::cout<<"elastic"<<std::endl;
  

  // Extract the components of the deformation field
  // -----------------------------------------------

  typename ComponentImageType::Pointer Fx = NULL;
  typename ComponentImageType::Pointer Fy = NULL;
  typename ComponentImageType::Pointer Fz = NULL;
  {
    typename VectorSplitterType::Pointer xCompExtractor = VectorSplitterType::New();
    xCompExtractor->SetInput( this->GetInput() );
    xCompExtractor->SetIndex(0);
    xCompExtractor->Update();
    Fx = xCompExtractor->GetOutput();

    typename VectorSplitterType::Pointer yCompExtractor = VectorSplitterType::New();
    yCompExtractor->SetInput( this->GetInput() );
    yCompExtractor->SetIndex(1);
    yCompExtractor->Update();
    Fy = yCompExtractor->GetOutput();

    typename VectorSplitterType::Pointer zCompExtractor = VectorSplitterType::New();
    zCompExtractor->SetInput( this->GetInput() );
    zCompExtractor->SetIndex(2);
    zCompExtractor->Update();
    Fz = zCompExtractor->GetOutput();
  }

  // Initial Gaussian smoothing
  // --------------------------

  if ( m_Verbose )
    std::cerr << "+ Gaussian smoothing, sigma=" << m_Sigma << "\n";

  typename ComponentImageType::Pointer Sx = NULL;
  typename ComponentImageType::Pointer Sy = NULL;
  typename ComponentImageType::Pointer Sz = NULL;

  {
    typename IIRSmootherType::Pointer smoothX = IIRSmootherType::New();
    smoothX->SetInput( Fx );
    smoothX->SetSigma( m_Sigma );
    smoothX->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
    smoothX->Update();
    Sx = smoothX->GetOutput();

    typename IIRSmootherType::Pointer smoothY = IIRSmootherType::New();
    smoothY->SetInput( Fy );
    smoothY->SetSigma( m_Sigma );
    smoothY->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
    smoothY->Update();
    Sy = smoothY->GetOutput();

    typename IIRSmootherType::Pointer smoothZ = IIRSmootherType::New();
    smoothZ->SetInput( Fz );
    smoothZ->SetSigma( m_Sigma );
    smoothZ->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
    smoothZ->Update();
    Sz = smoothZ->GetOutput();
  }

  // Elastic smoothing
  // -----------------

  if ( m_Kappa != 0.0 )
  {
    if ( m_Verbose )
      std::cerr << "+ Elastic smoothing, kappa=" << m_Kappa << "\n";

    typename ComponentImageType::SpacingType spacing = this->GetInput()->GetSpacing();

    // Caches used during the computation of the elastic smoothing
    typename ComponentImageType::Pointer Hx = ComponentImageType::New();
    Hx->CopyInformation( Fx );
    Hx->SetRegions( Fx->GetRequestedRegion() );
    Hx->Allocate();
    Hx->FillBuffer( 0.0 );

    typename ComponentImageType::Pointer Hy = ComponentImageType::New();
    Hy->CopyInformation( Fy );
    Hy->SetRegions( Fy->GetRequestedRegion() );
    Hy->Allocate();
    Hy->FillBuffer( 0.0 );

    typename ComponentImageType::Pointer Hz = ComponentImageType::New();
    Hz->CopyInformation( Fz );
    Hz->SetRegions( Fz->GetRequestedRegion() );
    Hz->Allocate();
    Hz->FillBuffer( 0.0 );

    unsigned int x,y,z;
    for ( unsigned int i = 0; i < 3; i++ )   // Rows of the Hessian
    {
      for ( unsigned int j = 0; j < 3; j++ ) // Columns of the Hessian
      {
        // Get the order of the derivation of each direction according to the
        // matrix position
        x = 0; y = 0; z = 0;

        switch(i)
        {
        case 0: x++; break;
        case 1: y++; break;
        case 2: z++; break;
        }

        switch(j)
        {
        case 0: x++; break;
        case 1: y++; break;
        case 2: z++; break;
        }

        // Initialize IIR filters
        typename IIRFilterType::Pointer xIIRFilter = IIRFilterType::New();
        xIIRFilter->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
        xIIRFilter->SetDirection( 0 );
        xIIRFilter->SetSigma( m_Sigma );
        switch(j)
        {
        case 0: xIIRFilter->SetInput( Fx ); break;
        case 1: xIIRFilter->SetInput( Fy ); break;
        case 2: xIIRFilter->SetInput( Fz ); break;
        }
        xIIRFilter->SetOrder( static_cast< typename IIRFilterType::OrderEnumType >(x) );

        typename IIRFilterType::Pointer yIIRFilter = IIRFilterType::New();
        yIIRFilter->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
        yIIRFilter->SetDirection( 1 );
        yIIRFilter->SetSigma( m_Sigma );
        yIIRFilter->SetInput( xIIRFilter->GetOutput() );
        yIIRFilter->SetOrder( static_cast< typename IIRFilterType::OrderEnumType >(y) );

        typename IIRFilterType::Pointer zIIRFilter = IIRFilterType::New();
        zIIRFilter->SetNormalizeAcrossScale( m_NormalizeAcrossScale );
        zIIRFilter->SetDirection( 2 );
        zIIRFilter->SetSigma( m_Sigma );
        zIIRFilter->SetInput( yIIRFilter->GetOutput() );
        zIIRFilter->SetOrder( static_cast< typename IIRFilterType::OrderEnumType >(z) );

        // Compute spacing coefficient
        double coef = 1.0;
        if ( x != 0 )
          coef *= pow( spacing[0], (double)x );
        if ( y != 0 )
          coef *= pow( spacing[1], (double)y );
        if ( z != 0 )
          coef *= pow( spacing[2], (double)z );

        // Multiply by the spacing coefficient
        typename DividerType::Pointer divider = DividerType::New();
        divider->SetInput( zIIRFilter->GetOutput() );
        divider->SetConstant( coef );

        // Adding (in-place)
        typename AdderType::Pointer adder = AdderType::New();
        adder->SetInput1( divider->GetOutput() );
        adder->InPlaceOn();
        switch(i)
        {
        case 0:
          adder->SetInput2( Hx );
          adder->Update();
          Hx = adder->GetOutput();
          break;
        case 1:
          adder->SetInput2( Hy );
          adder->Update();
          Hy = adder->GetOutput();
          break;
        case 2:
          adder->SetInput2( Hz );
          adder->Update();
          Hz = adder->GetOutput();
          break;
        }

      } // End loop on the columns j
    } // End loop on the lines i

    // Mask the elastic component
    if ( m_MaskElasticSmoothing && m_MaskImage )
    {
      if ( m_Verbose )
        std::cerr << "+ Masking\n";
      {
        typename MaskerType::Pointer masker = MaskerType::New();
        masker->SetInput1( Hx );
        masker->SetInput2( m_MaskImage );
        masker->Update();
        Hx = masker->GetOutput();
      }

      {
        typename MaskerType::Pointer masker = MaskerType::New();
        masker->SetInput1( Hy );
        masker->SetInput2( m_MaskImage );
        masker->Update();
        Hy = masker->GetOutput();
      }

      {
        typename MaskerType::Pointer masker = MaskerType::New();
        masker->SetInput1( Hz );
        masker->SetInput2( m_MaskImage );
        masker->Update();
        Hz = masker->GetOutput();
      }
    }

    // Add the resulting elastic component region to the Gaussian filter
    {
      typename MultiplierType::Pointer multiplier = MultiplierType::New();
      multiplier->SetInput( Hx );
      multiplier->SetConstant( m_Sigma * m_Sigma * m_Kappa / (1.0 + m_Kappa) );

      typename AdderType::Pointer adder = AdderType::New();
      adder->SetInput1( Sx );
      adder->SetInput2( multiplier->GetOutput() );
      adder->Update();
      Sx = adder->GetOutput();
    }

    {
      typename MultiplierType::Pointer multiplier = MultiplierType::New();
      multiplier->SetInput( Hy );
      multiplier->SetConstant( m_Sigma * m_Sigma * m_Kappa / (1.0 + m_Kappa) );

      typename AdderType::Pointer adder = AdderType::New();
      adder->SetInput1( Sy );
      adder->SetInput2( multiplier->GetOutput() );
      adder->Update();
      Sy = adder->GetOutput();
    }

    {
      typename MultiplierType::Pointer multiplier = MultiplierType::New();
      multiplier->SetInput( Hz );
      multiplier->SetConstant( m_Sigma * m_Sigma * m_Kappa / (1.0 + m_Kappa) );

      typename AdderType::Pointer adder = AdderType::New();
      adder->SetInput1( Sz );
      adder->SetInput2( multiplier->GetOutput() );
      adder->Update();
      Sz = adder->GetOutput();
    }

  } // End elastic smoothing

  // Rebuild the vector field
  typename VectorComposerType::Pointer ComposeVectorFilter = VectorComposerType::New();
  ComposeVectorFilter->SetInput1( Sx );
  ComposeVectorFilter->SetInput2( Sy );
  ComposeVectorFilter->SetInput3( Sz );
  ComposeVectorFilter->GraftOutput( this->GetOutput() );
  ComposeVectorFilter->Update();

  this->GraftOutput( ComposeVectorFilter->GetOutput() );
  this->GetOutput()->Modified();
}



// -----------------------------------------------------------------------------
// PrintSelf method

template <class TImageType, class TMaskImageType >
void
ElasticVectorSmoothingFilter< TImageType, TMaskImageType >::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Sigma: "
      << m_Sigma << std::endl;

  os << indent << "Kappa: "
      << m_Kappa << std::endl;

  os << indent << "Mask elastic smoothing: "
      << m_MaskElasticSmoothing << std::endl;

  os << indent << "Normalize across scale: "
      << m_NormalizeAcrossScale << std::endl;

  os << indent << "Verbose: "
      << m_Verbose << std::endl;
}



} // end namespace itk

#endif

