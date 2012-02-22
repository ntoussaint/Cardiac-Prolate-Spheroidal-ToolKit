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
#include <itkResampleImage3Command.h>

#include "itkLinearInterpolateImageFunction.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include "GetPot.h"

namespace itk
{

  ResampleImage3Command::ResampleImage3Command()
  {
    m_ShortDescription = "Resample a 3D image to a reference field-of-view";
    m_LongDescription = m_ShortDescription;
    m_LongDescription = "\n\nUsage:\n";
    m_LongDescription +="-i  [input image]\n";    
    m_LongDescription +="-r  [reference field-of-view image ]\n";    
    m_LongDescription +="-o  [output image]\n";    
  }

  ResampleImage3Command::~ResampleImage3Command()
  {}

  int ResampleImage3Command::Execute (int narg, const char* arg[])
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
  
    const char* File = cl.follow("NoFile",2,"-I","-i");
    const char* ref        = cl.follow("NoFile",2,"-r","-R");
    const char* outFile    = cl.follow("NoFile",2,"-O","-o");
    
    typedef double                                ScalarType;  
    typedef itk::Image<ScalarType, 3> ImageType;
    typedef itk::ImageFileReader <ImageType> ImageReaderType;
    typedef itk::ImageFileWriter <ImageType> ImageWriterType;
    typedef ImageType::SizeType    SizeType;
    typedef ImageType::SpacingType SpacingType;
    typedef ImageType::PointType   PointType;
    typedef itk::ImageRegionIterator<ImageType> ImageIteratorType;
    typedef itk::LinearInterpolateImageFunction<ImageType, double>  InterpolatorType;
  
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(File);
    std::cout << "Reading: " << File << std::endl;
    try
    {
      reader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      exit(-1);
    }
    ImageType::Pointer image = reader->GetOutput();
  
    ImageReaderType::Pointer refreader = ImageReaderType::New();
    refreader->SetFileName(ref);
    std::cout << "Reading: " << ref << std::endl;
    try
    {
      refreader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      exit(-1);
    }
    ImageType::Pointer reference = refreader->GetOutput();

  
    std::cout<<" allocating image."<<std::endl;
    ImageType::Pointer output = ImageType::New();
    output->SetRegions (reference->GetLargestPossibleRegion());
    output->SetOrigin(reference->GetOrigin());
    output->SetSpacing(reference->GetSpacing());  
    output->SetDirection(reference->GetDirection());
    output->Allocate();

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetInputImage(image);

    itk::ImageRegionIterator<ImageType> itIn(image, image->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType> itOut(output, output->GetLargestPossibleRegion());
    ImageType::PointType x;
    itk::ContinuousIndex<double, 3> index;
    itIn.GoToBegin();
    itOut.GoToBegin();
    unsigned int counter = 0;
  
    std::cout<<" iterate."<<std::endl;
    while( !itOut.IsAtEnd() )
    {
      output->TransformIndexToPhysicalPoint (itOut.GetIndex(), x);
      ScalarType value = static_cast<ScalarType>(0.0);
      bool isinside = image->TransformPhysicalPointToContinuousIndex (x, index);
      if (isinside)
	value = interpolator->EvaluateAtContinuousIndex (index);
      else
	counter++;
    
      itOut.Set (value);
      ++itOut;
    }
  
    ImageWriterType::Pointer writer = ImageWriterType::New();
  
    writer->SetFileName(outFile);
    writer->SetInput( output );
  
    std::cout << "Writing: " << outFile << std::flush;
    try
    {
      writer->Update();
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
