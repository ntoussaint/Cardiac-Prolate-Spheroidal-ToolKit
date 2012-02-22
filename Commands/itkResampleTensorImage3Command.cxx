/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkResampleTensorImage3.cxx 1 2010-05-21 14:00:33Z nt08 $
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

/**
 * 
 * \file itkResampleTensorImage3.cxx
 *
 * \brief Yet another version of resamping tensors
 * 
 * This function takes a tensor field as an input in ITK format. It can take
 * one or several other inputs.
 * \param matrix an rigid or affine transform to transform tensors
 * \param reference The reference image to use as the output space
 * \param sx The size of the output in pixels in the three directions
 * Then the origin directions and spacing of either the ref or the input are used. 
 *
 * 
 * \author Nicolas Toussaint, King's College London, nicolas.toussaint@kcl.ac.uk
 *
 * \ingroup TensorProcessing
 * 
 * \url http://www.kcl.ac.uk/schools/medicine/research/imaging/
 *  
 */
#include <itkResampleTensorImage3Command.h>

#include "itkLogTensorImageFilter.h"
#include "itkExpTensorImageFilter.h"
#include "itkResampleTensorImageFilter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkAffineTensorTransform.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include "itkTensorLinearInterpolateImageFunction.h"
#include "itkTensorImageIO.h"
#include "GetPot.h"


namespace itk
{

  ResampleTensorImage3Command::ResampleTensorImage3Command()
  {
    m_ShortDescription = "Resample a 3D tensor image to a reference field-of-view";
    m_LongDescription = m_ShortDescription;
    m_LongDescription = "\n\nUsage:\n";
    m_LongDescription +="-i  [input image]\n";    
    m_LongDescription +="-r  [reference field-of-view image ]\n";    
    m_LongDescription +="-o  [output image]\n";    
  }

  ResampleTensorImage3Command::~ResampleTensorImage3Command()
  {}

  int ResampleTensorImage3Command::Execute (int narg, const char* arg[])
  {

    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const bool IsInputPresent = cl.search(2,"-I","-i");
    const bool IsOutputPresent = cl.search(2,"-O","-o");

    if( !IsInputPresent || !IsOutputPresent )
    {
      std::cerr << "Error: Input and (or) output not set." << std::endl;
      exit (-1);
    }
  
    const char* tensorFile = cl.follow("NoFile",2,"-I","-i");
    const char* ref        = cl.follow("NoFile",2,"-r","-R");
    const int sx           = cl.follow(-1,2,"-sx","-SX");
    const int sy           = cl.follow(-1,2,"-sy","-SX");
    const int sz           = cl.follow(-1,2,"-sz","-SX");
    const char* outFile    = cl.follow("NoFile",2,"-O","-o");
    const char* mat        = cl.follow("NoFile",2,"-m","-M");
    
    typedef double                                ScalarType;  
    typedef itk::TensorImageIO<ScalarType, 3, 3>  IOType;
    typedef IOType::TensorImageType               TensorImageType;
    typedef itk::ResampleTensorImageFilter<TensorImageType,TensorImageType> FilterType;
    typedef TensorImageType::SizeType    SizeType;
    typedef TensorImageType::SpacingType SpacingType;
    typedef TensorImageType::PointType   PointType;

    typedef itk::MatrixOffsetTransformBase<ScalarType, 3, 3> TransformType;
    typedef itk::AffineTensorTransform< double, 3 >  TensorTransformType;
    typedef itk::TransformFileReader TransformReaderType;
    // register linear transform into factory for I/O purposes
    itk::TransformFactory<TransformType>::RegisterTransform ();
    // Software Guide : EndCodeSnippet

    IOType::Pointer reader = IOType::New();
    reader->SetFileName(tensorFile);

    std::cout << "Reading: " << tensorFile << std::endl;
    try
    {
      reader->Read();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      exit(-1);
    }
  
    TensorImageType::Pointer tensors = reader->GetOutput();

    std::cout<<tensors<<std::endl;
  

    std::cout<<"entering the transform file"<<std::endl;
  
    TensorTransformType::Pointer tensortransform = TensorTransformType::New();

    if (strcmp (mat, "NoFile"))
    {
    
      TransformReaderType::Pointer transformReader
	= TransformReaderType::New();    
      transformReader->SetFileName(  mat );
    
      // Update the reader
      try
      {
	transformReader->Update();
      }
      catch( itk::ExceptionObject& err )
      {
	std::cout << "Could not read the input transform." << std::endl;
	std::cout << err << std::endl;
	exit( EXIT_FAILURE );
      }
    
      typedef TransformReaderType::TransformType BaseTransformType;
      BaseTransformType* baseTrsf(0);
      const TransformReaderType::TransformListType* trsflistptr
	= transformReader->GetTransformList();
      if ( trsflistptr->empty() )
      {
	std::cout << "Could not read the input transform." << std::endl;
	exit( EXIT_FAILURE );
      }
      else if (trsflistptr->size()>1 )
      {
	std::cout << "The input transform file contains more than one transform, we use the first one." << std::endl;
      }
    
      baseTrsf = trsflistptr->front();
      if ( !baseTrsf )
      {
	std::cout << "Could not read the input transform." << std::endl;
	exit( EXIT_FAILURE );
      }
    
      TransformType::Pointer transform = dynamic_cast<TransformType*>(baseTrsf);
      if ( !transform )
      {
	std::cout << "Could not cast input transform to a usable transform for the factory." << std::endl;
	exit( EXIT_FAILURE );
      }
    
      std::cout << transform << std::endl;

      tensortransform->SetCenter (transform->GetCenter());
      tensortransform->SetTranslation (transform->GetTranslation());
      tensortransform->SetMatrix (transform->GetMatrix());

    }
  
    FilterType::Pointer filter = FilterType::New();

    typedef itk::TensorLinearInterpolateImageFunction<TensorImageType, double>  InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
  
    filter->SetTensorInterpolator( interpolator );
    filter->SetInput ( tensors );
    if (strcmp (mat, "NoFile"))
      filter->SetTensorTransform( tensortransform );
  
    TensorImageType::SpacingType spacing;
    TensorImageType::PointType origin;
    TensorImageType::SizeType size;
    TensorImageType::DirectionType direction;

    if (strcmp (ref, "NoFile"))
    {
    
      typedef itk::Image<short, 3> ImageType;
      ImageType::Pointer reference;
      itk::ImageFileReader< ImageType >::Pointer io2 = itk::ImageFileReader< ImageType >::New();
      io2->SetFileName( ref );
      try
      {
	io2->Update();
      } catch (itk::ExceptionObject &e)
      {
	std::cerr << e;
	return -1;
      }

      reference = io2->GetOutput();

      spacing = reference->GetSpacing();
      origin  = reference->GetOrigin();
      size =   reference->GetLargestPossibleRegion().GetSize();
      direction = reference->GetDirection();
    }
    else if ( (sx != -1) && (sy != -1) && (sz != -1) )
    {
      spacing = tensors->GetSpacing();
      origin  = tensors->GetOrigin();
      size =   tensors->GetLargestPossibleRegion().GetSize();
      direction = tensors->GetDirection();
    
      spacing[0] /= (double)(sx) / (double)(size[0]);
      spacing[1] /= (double)(sy) / (double)(size[1]);
      spacing[2] /= (double)(sz) / (double)(size[2]);

      size[0] = sx;
      size[1] = sy;
      size[2] = sz;    
    }
    else
    {
      spacing = tensors->GetSpacing();
      origin  = tensors->GetOrigin();
      size =   tensors->GetLargestPossibleRegion().GetSize();
      direction = tensors->GetDirection();
    }
  
    filter->SetOutputOrigin( origin );
    filter->SetOutputSpacing( spacing );
    filter->SetSize( size );
    filter->SetOutputDirection( direction );

    try
    {
      filter->Update();
    }
    catch (itk::ExceptionObject &e )
    {
      std::cerr << e;
      exit (-1);
    }

    tensors = filter->GetOutput();  
    filter->GetOutput()->Update();

    IOType::Pointer writer = IOType::New();
  
    writer->SetFileName(outFile);
    writer->SetInput( filter->GetOutput() );
  
    std::cout << "Writing: " << outFile << std::flush;
    try
    {
      writer->Write();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      exit(-1);
    }
    std::cout << " Done." << std::endl;
  
    return EXIT_SUCCESS;

  
  }
}
