/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkApplyTransformToImage.cxx 1 2010-05-21 14:00:33Z nt08 $
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
 * \file itkApplyTransformToImage.cxx
 *
 * \brief  Apply a simple transform to an image, in ITK format. 
 * 
 * This tool takes an image and a transform in ITK format as inputs.
 * The transform is applied to the "header" of the image, meaning
 * it will only transform the origin and 3x3 direction of the input
 * image. Therefore, the input transform has to be a derivative of
 * a MatrixOffsetTransformBase, i.e. an affine transform
 *
 * \note This tool does *NOT* check if the input transform matrix
 * is rigid. If it is not the case the output image could have direction
 * vectors that are *NOT* orthogonal to eachother. 
 * 
 * \author Nicolas Toussaint, King's College London, nicolas.toussaint@kcl.ac.uk
 * 
 * \url http://www.kcl.ac.uk/schools/medicine/research/imaging/
 *  
 */
#include <itkApplyTransformToImageCommand.h>

#include <itkImage.h>
#include <cstdlib>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkTransformFileReader.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkResampleImageFilter.h>
#include <itkTransformFactory.h>
#include "itkTranslationTransform.h"

#include "GetPot.h"


namespace itk
{

  template<typename TranslationTransformType, typename UnderImageType, typename ImageType>
  void ApplyResampledPixelContainer(const char* file, typename TranslationTransformType::Pointer transform, typename ImageType::Pointer image)
  {
    typedef ImageFileReader<UnderImageType> ReaderType;
    typedef ResampleImageFilter< UnderImageType, UnderImageType > ResampleFilterType;
    
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName (file);
    reader->Update();
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput (reader->GetOutput());
    resampler->SetTransform (transform);
    resampler->SetDefaultPixelValue( 0 );
    resampler->SetOutputParametersFromImage (reader->GetOutput());
    resampler->Update();
    
    image->SetPixelContainer (resampler->GetOutput()->GetPixelContainer());
  }

  ApplyTransformToImageCommand::ApplyTransformToImageCommand()
  {
    m_ShortDescription = "Apply a rigid transformation to an image without resampling";
    m_LongDescription = m_ShortDescription;
    m_LongDescription = "\n\nUsage:\n";
    m_LongDescription +="-i [input image]\n";
    m_LongDescription +="-t [transform]\n";
    m_LongDescription +="-o [Output file name]\n";
  }

  ApplyTransformToImageCommand::~ApplyTransformToImageCommand()
  {}

  int ApplyTransformToImageCommand::Execute (int narg, const char* arg[])
  {


    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const bool IsInputPresent    = cl.search(2,"-i","-I");
    const bool IsOutputPresent   = cl.search(2,"-o","-O");
  
    if(!IsInputPresent || !IsOutputPresent )
    {
      std::cerr << "Error: Input and (or) output not set." << std::endl;
      exit (EXIT_FAILURE);
    }


    const char* fileOut  = cl.follow("NoFile",2,"-o","-O");
    const char* fileIn  = cl.follow("NoFile",2,"-i","-I");
    const char* fileTr  = cl.follow("NoFile",2,"-t","-T");

    
    static const unsigned int Dimension = 3;
    static const unsigned int UnderDimension = 2;

    typedef double                                       ScalarType;  
    typedef itk::Image<ScalarType, Dimension>                    ImageType;
    typedef ImageType::DirectionType                     DirectionType;

    typedef itk::MatrixOffsetTransformBase<double, Dimension, Dimension> LinearTransformType;
    typedef itk::TranslationTransform<double, UnderDimension>         UnderTranslationTransformType;
    typedef itk::TranslationTransform<double, Dimension>         TranslationTransformType;
    typedef itk::Transform<double, Dimension, Dimension>                 TransformType;
    typedef itk::TransformFileReader                     TransformReaderType;
    typedef ImageType::PointType                         PointType;
    typedef itk::ImageFileReader <ImageType>             ImageReaderType;  
    typedef itk::ImageFileWriter <ImageType>             ImageWriterType;
    typedef itk::Image<ScalarType, UnderDimension> UnderImageType;
    typedef itk::ImageFileReader <UnderImageType>  UnderImageReaderType; 
    typedef itk::ResampleImageFilter< UnderImageType, UnderImageType > ResampleFilterType;
    typedef UnderImageType::PixelContainer PixelContainerType;

    itk::TransformFactory< LinearTransformType >::RegisterTransform ();
    itk::TransformFactory< UnderTranslationTransformType >::RegisterTransform ();
    itk::TransformFactory< TranslationTransformType >::RegisterTransform ();
  
    std::cout<<"reading input file "<<fileIn<<std::endl;  
    ImageReaderType::Pointer myReader = ImageReaderType::New();
    myReader->SetFileName (fileIn);
    try
    {
      myReader->Update();
    } catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return EXIT_FAILURE;
    }  
    ImageType::Pointer myImage = myReader->GetOutput();
    std::cout<<"done."<<std::endl;

  
    TransformReaderType::Pointer reader = TransformReaderType::New();
    reader->SetFileName ( fileTr );
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return EXIT_FAILURE;
    }
    
    LinearTransformType::Pointer      transform1 = dynamic_cast<LinearTransformType*>( reader->GetTransformList()->front().GetPointer() );
    TranslationTransformType::Pointer transform2 = dynamic_cast<TranslationTransformType*>( reader->GetTransformList()->front().GetPointer() );
    UnderTranslationTransformType::Pointer transform3 = dynamic_cast<UnderTranslationTransformType*>( reader->GetTransformList()->front().GetPointer() );
    
    if (!transform1 && !transform2 && !transform3)
    {
      std::cerr << "The transformation written in "<<fileTr<<" does not derive from a MatrixOffsetTransformBase, "
		<<"which is mandatory in this executable to be able to only change the header of the file"<<std::endl;
      return EXIT_FAILURE;
    }
    
    PointType origin = myImage->GetOrigin();
    PointType neworigin = origin;
    if (transform1)
      neworigin = transform1->TransformPoint (origin);
    else if (transform2)
      neworigin = transform2->TransformPoint (origin);
    
    DirectionType direction = myImage->GetDirection();
    DirectionType newdirection = direction;
    if (transform1)
      newdirection = transform1->GetMatrix() * direction;
    else if (transform2)
      newdirection = direction;
  
    myImage->SetOrigin(neworigin);
    myImage->SetDirection(newdirection);

    if (transform3)
    {
      std::cout<<"detected a 2D transform : "<<std::endl;
      transform3->Print (std::cout);
      ApplyResampledPixelContainer<UnderTranslationTransformType,UnderImageType, ImageType>(fileIn, transform3, myImage);
    }
  
    std::cout<<"writing output file "<<fileOut<<std::endl;  
    ImageWriterType::Pointer myWriter = ImageWriterType::New();
    myWriter->SetFileName (fileOut);
    myWriter->SetInput (myImage);
  
    try
    {
      myWriter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return EXIT_FAILURE;
    }

  
    std::cout << " Done." << std::endl;  
  
    return EXIT_SUCCESS;

  }
}
