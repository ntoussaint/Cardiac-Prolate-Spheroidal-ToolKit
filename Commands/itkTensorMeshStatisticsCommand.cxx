/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkTensorMeshToImage.cxx 1 2010-05-21 14:00:33Z nt08 $
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
#include <itkTensorMeshStatisticsCommand.h>

#include <itkTensorMeshStatistics.h>

#include <itkTensorImageIO.h>

#include <itkTensorMeshIO.h>

#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkGradientTensorImageFilter.h>
#include <itkEigenVectorToTensorImageFilter.h>
#include <itkLogTensorImageFilter.h>

#include "GetPot.h"

namespace itk
{

  TensorMeshStatisticsCommand::TensorMeshStatisticsCommand()
  {
    m_ShortDescription = "Compute statistics in prolate sph. coordinates over a set of tensors";
    m_LongDescription += m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-i    [input  tensor image]\n";
    m_LongDescription += "-pr   [prolate transform used]\n";
    m_LongDescription += "-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription += "-f2   [backward displacement field (default : backward.mha)]\n";  
    m_LongDescription += "-s    [statistic output type]\n";  
    m_LongDescription += "      [1: Norm of the CovarianceMatrix : variability metric]\n";
    m_LongDescription += "      [2: FA (Fractional Anisotropy) of the CovarianceMatrix : variability anisotropy metric]\n";
    m_LongDescription += "      [3: ]\n";
    m_LongDescription += "-o    [output image]\n";
    
  }

  TensorMeshStatisticsCommand::~TensorMeshStatisticsCommand()
  {}

  int TensorMeshStatisticsCommand::Execute (int narg, const char* arg[])
  {
    
    typedef double                                      ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>        TensorIOType;
    typedef itk::TensorMeshIO<ScalarType, 3, 3>         TensorMeshIOType;
    typedef itk::TensorMeshStatistics<ScalarType, 3>    FilterType;
    typedef FilterType::TensorType                      TensorType;
    typedef FilterType::TransformType                   TransformType;
    typedef FilterType::InputImageType                  TensorImageType;
    typedef FilterType::OutputImageType                 ImageType;
    typedef FilterType::DisplacementFieldType           DisplacementFieldType;
    typedef FilterType::MeshType                        MeshType;
    typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFileReaderType;
    typedef itk::ImageFileWriter<ImageType>             ImageFileWriterType;
    typedef itk::GradientTensorImageFilter<TensorImageType, ScalarType> GradientFilterType;
    typedef GradientFilterType::OutputImageType EigenVectorImageType;
    typedef itk::EigenVectorToTensorImageFilter<EigenVectorImageType,TensorImageType> VectorToTensorFilterType;
    typedef itk::LogTensorImageFilter<TensorImageType, TensorImageType> LogFilterType;
 
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    const char* inputfile                    = cl.follow("input.vtk",2,"-i","-I");
    const char* prolatefile                  = cl.follow("prolate.pr",2,"-pr","-PR");
    const char* displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char* outputfile                   = cl.follow("output.mha",2,"-o","-O");
    const unsigned int statistictype         = cl.follow(1,2,"-s","-S");
  
    std::cout << "inputfile: \t\t\t" << inputfile << std::endl;
    std::cout << "prolatefile: " << prolatefile << std::endl;
    std::cout << "displacementfieldfile: " << displacementfieldfile << std::endl;
    std::cout << "inversedisplacementfieldfile: " << inversedisplacementfieldfile << std::endl;
    std::cout << "outputfile: \t\t\t" <<outputfile << std::endl;
    std::cout << "statistictype:\t\t " <<statistictype<<" [";
    switch(statistictype)
    {
	case 1:
	  std::cout << "Norm of the CovarianceMatrix]"<<std::endl;
	  break;
	case 2:
	  std::cout << "FA (Fractional Anisotropy) of the CovarianceMatrix]"<<std::endl;
	  break;
	default:
	  std::cout << "Unknown type]"<<std::endl;
	  std::cout << std::endl << this->GetLongDescription() << std::endl;
	  return -1;
	  break;
    }
  
    std::cout << std::flush;
  
    std::cout << "Reading input tensor image: " << inputfile <<"... "<< std::flush;
    TensorIOType::Pointer reader = TensorIOType::New();
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

    // LogFilterType::Pointer logger = LogFilterType::New();
    // logger->SetInput (reader->GetOutput());
  
    // try
    // {
    //   logger->Update();
    // }
    // catch(itk::ExceptionObject &e)
    // {
    //   std::cerr << e << std::endl;
    //   std::exit(EXIT_FAILURE);
    // }
    // std::cout << " Done." << std::endl;
  
    // GradientFilterType::Pointer gradientfilter = GradientFilterType::New();
    // gradientfilter->SetInput (logger->GetOutput());
  
    // try
    // {
    //   gradientfilter->Update();
    // }
    // catch(itk::ExceptionObject &e)
    // {
    //   std::cerr << e << std::endl;
    //   std::exit(EXIT_FAILURE);
    // }
    // std::cout << " Done." << std::endl;

    // VectorToTensorFilterType::Pointer vector2tensor = VectorToTensorFilterType::New();
    // vector2tensor->SetInput (gradientfilter->GetOutput());
  
    // try
    // {
    //   vector2tensor->Update();
    // }
    // catch(itk::ExceptionObject &e)
    // {
    //   std::cerr << e << std::endl;
    //   std::exit(EXIT_FAILURE);
    // }
    // std::cout << " Done." << std::endl;

    // reader->SetFileName("test.mha");
    // reader->SetInput (vector2tensor->GetOutput());
    // try
    // {
    //   reader->Write();
    // }
    // catch(itk::ExceptionObject &e)
    // {
    //   std::cerr << e << std::endl;
    //   std::exit(EXIT_FAILURE);
    // }
    // std::cout << " Done." << std::endl;


  
    std::cout << "Reading forward field: " << displacementfieldfile <<"... "<<std::flush;
    DisplacementFileReaderType::Pointer displacementreader = DisplacementFileReaderType::New();
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
    DisplacementFieldType::Pointer displacementfield = displacementreader->GetOutput();
    std::cout << " Done." << std::endl;

    std::cout << "Reading backward field: " << inversedisplacementfieldfile <<"... "<<std::flush;
    DisplacementFileReaderType::Pointer displacementreader2 = DisplacementFileReaderType::New();
    displacementreader2->SetFileName(inversedisplacementfieldfile);
    try
    {
      displacementreader2->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    DisplacementFieldType::Pointer inversedisplacementfield = displacementreader2->GetOutput();
    std::cout << " Done." << std::endl;

    std::cout<<"reading transform " << prolatefile <<"... "<<std::flush;
    itk::TransformFactory<TransformType>::RegisterTransform ();  
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
    TransformType::Pointer transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    std::cout << " Done." << std::endl;

    FilterType::Pointer statistics = FilterType::New();
    switch(statistictype)
    {
	case 2:
	  statistics->SetStatisticsOutputType (FilterType::CovarianceMatrixFA);
	  break;
	case 1:
	default:
	  statistics->SetStatisticsOutputType (FilterType::CovarianceMatrixNorm);
	  break;
    }
  
    statistics->SetInput (reader->GetOutput());
    statistics->SetTransform (transform);
    statistics->SetDisplacementField (displacementfield);
    statistics->SetInverseDisplacementField (inversedisplacementfield);

    statistics->GetLimiter()->CanineDivisionsOn();
    statistics->GetLimiter()->SetAHASegmentationType (FilterType::AHALimiterType::AHA_17_ZONES);

    MeshType::Pointer outputmesh = statistics->GetMeshOutput();
    TensorMeshIOType::Pointer meshwriter = TensorMeshIOType::New();
    meshwriter->SetInput (outputmesh);
  
    try
    {
      statistics->SetCoordinatesType (FilterType::ProlateSpheroidal);
      statistics->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
  
    ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
    writer->SetFileName (outputfile);
    writer->SetInput (statistics->GetOutput());
    meshwriter->SetFileName ("outputmesh.vtk");
  
    try
    {
      writer->Update();
      meshwriter->Write();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
  
  }

}

  
