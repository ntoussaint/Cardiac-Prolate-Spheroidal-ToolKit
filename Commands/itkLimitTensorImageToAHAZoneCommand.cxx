/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkExtrapolateTensorField.cxx 1 2010-05-21 14:00:33Z nt08 $
  Language:  C++
  Author:    $Author: nt08 $
  Date:      $Date: 2010-05-21 14:00:33 +0000 (Fri, 21 May 2010) $
  Version:   $Revision: 1 $

  Copyright (c) 2010 King's College London - Division of Imaging Sciences. All rights reserved.
  See Copyright.txt for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

  ============================================================================*/
#include <itkLimitTensorImageToAHAZoneCommand.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkLimitToAHAZoneImageFilter.h>
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

#include "itkTensorImageIO.h"
#include "itkTensorMeshIO.h"

#include "GetPot.h"

namespace itk
{

  LimitTensorImageToAHAZoneCommand::LimitTensorImageToAHAZoneCommand()
  {
    m_ShortDescription = "Crop a a tensor image/mesh to an AHA zone in the Prolate Spheroidal sense";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input tensor image/mesh]\n";
    m_LongDescription +="-f    [inverse displacement field (default : backward.mha)]\n";
    m_LongDescription +="-pr   [prolate transform used (default : prolate.pr)]\n";
    m_LongDescription +="-z    [AHA zone to take into account (1~17) (default: 0 : ALL ZONES)]\n";
    m_LongDescription +="-o    [output image zone file]\n";
  }

  LimitTensorImageToAHAZoneCommand::~LimitTensorImageToAHAZoneCommand()
  {}

  int LimitTensorImageToAHAZoneCommand::Execute (int narg, const char* arg[])
  {
  
    typedef double                               ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3> TensorImageIOType;
    typedef itk::TensorMeshIO<ScalarType, 3, 3>  TensorMeshIOType;
    typedef TensorImageIOType::TensorImageType        ImageType;  
    typedef itk::Vector<float, 3>                DisplacementType;
    typedef itk::Image<DisplacementType, 3>      DisplacementFieldType;
    typedef itk::ProlateSpheroidalTransform<ScalarType> TransformType;
    typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFileReaderType;
    typedef itk::LimitToAHAZoneImageFilter<ImageType> AHALimiterType;    

    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    const char*  inputfile                    = cl.follow("input.vtk",2,"-i","-I");
    const char*  prolatefile                  = cl.follow("prolate.pr",2,"-pr","-PR");
    const char*  displacementfieldfile        = cl.follow("forward.mha",2,"-f","-F");
    const char*  outputfile                   = cl.follow("output.csv",2,"-o","-O");
    const unsigned int zone                   = cl.follow(0, 2,"-z","-Z");

    std::cout << "Processing AHA segments: " << std::endl;
    std::cout << "inputfile: \t\t\t" << inputfile << std::endl;
    std::cout << "prolatefile: \t\t\t" << prolatefile << std::endl;
    std::cout << "displacementfieldfile: \t\t" << displacementfieldfile << std::endl;
    std::cout << "outputfile: \t\t\t" <<outputfile << std::endl;
    std::cout << "zone:      \t\t\t" <<zone << std::endl;
  
    std::cout << std::flush;

    std::cout << "Reading input image: " << inputfile << std::flush;  
    TensorImageIOType::Pointer reader = TensorImageIOType::New();
    reader->SetFileName(inputfile);
    try
    {
      reader->Read();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << " Done." << std::endl;
    ImageType::Pointer input = reader->GetOutput();
  
    std::cout<<"reading prolate transform"<<std::endl;
    TransformType::Pointer transform = TransformType::New();  
    itk::TransformFactory<TransformType>::RegisterTransform ();
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
    transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    std::cout << " Done." << std::endl;

    std::cout << "Reading forward field: " << displacementfieldfile << std::flush;  
    DisplacementFileReaderType::Pointer    displacementreader = DisplacementFileReaderType::New();
    displacementreader->SetFileName(displacementfieldfile);
    try
    {
      displacementreader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << " Done." << std::endl;
    DisplacementFieldType::Pointer displacementfield = displacementreader->GetOutput();

    std::cout<<"limiting the image "<<inputfile<<" to AHA zone: "<<zone<<std::endl;
    AHALimiterType::Pointer zonelimiter = AHALimiterType::New();
    zonelimiter->SetInput (input);
    zonelimiter->SetAHAZone (zone);
    zonelimiter->SetTransform (transform);
    zonelimiter->SetInverseDisplacementField (displacementfield);
    zonelimiter->CanineDivisionsOff();
    zonelimiter->SetAHASegmentationType (AHALimiterType::AHA_17_ZONES);
    // zonelimiter->SetAHASegmentationType (AHALimiterType::AHA_3_ZONES);
    zonelimiter->CalculateZones();

    try
    {
      zonelimiter->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout<<"Done."<<std::endl;
    ImageType::Pointer output = zonelimiter->GetOutput();

    std::cout<<"writing resuting image to : "<<outputfile<<std::endl;
    TensorImageIOType::Pointer writer = TensorImageIOType::New();
    writer->SetFileName (outputfile);
    writer->SetInput (output);
    try
    {
      writer->Write();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
  
  }

  
}
