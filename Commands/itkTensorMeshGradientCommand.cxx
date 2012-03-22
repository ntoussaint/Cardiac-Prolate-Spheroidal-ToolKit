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
#include <itkTensorMeshGradientCommand.h>

#include <itkTensorImageIO.h>
#include <itkTensorMeshIO.h>

#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "itkProlateSpheroidalGradientTensorMeshFilter.h"
#include "itkTensorImageToMeshFilter.h"
#include "itkTensorMeshToImageFilter.h"
#include "itkWarpTensorMeshFilter.h"
#include "itkProlateSpheroidalTransformTensorMeshFilter.h"

#include "GetPot.h"

namespace itk
{

  TensorMeshGradientCommand::TensorMeshGradientCommand()
  {
    m_ShortDescription = "Compute the gradient of a tensor field in Prolate Spheroidal sense";
    m_LongDescription += m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-i    [input  tensor image]\n";
    m_LongDescription += "-pr   [prolate transform used]\n";
    m_LongDescription += "-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription += "-f2   [backward displacement field (default : backward.mha)]\n";  
    m_LongDescription += "-o    [output image]\n";
    
  }

  TensorMeshGradientCommand::~TensorMeshGradientCommand()
  {}

  int TensorMeshGradientCommand::Execute (int narg, const char* arg[])
  {
    
    typedef double                                      ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>        TensorIOType;
    typedef itk::TensorMeshIO<ScalarType, 3, 3>         TensorMeshIOType;
    typedef TensorMeshIOType::TensorMeshType            MeshType;
    
    typedef itk::ProlateSpheroidalGradientTensorMeshFilter<MeshType> FilterType;
    typedef FilterType::TensorType                      TensorType;
    typedef FilterType::TransformType                   TransformType;
    typedef TensorIOType::TensorImageType               TensorImageType;
    typedef itk::Vector<ScalarType, 3>                  DisplacementType;
    typedef itk::Image<DisplacementType, 3>             DisplacementFieldType;
    typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFileReaderType;
    typedef itk::TensorImageToMeshFilter<TensorType, 3> ImageToMeshFilterType;

    typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>  WarperType;
    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>   TransformerType;    
 
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
  
    std::cout << "inputfile: \t\t\t" << inputfile << std::endl;
    std::cout << "prolatefile: " << prolatefile << std::endl;
    std::cout << "displacementfieldfile: " << displacementfieldfile << std::endl;
    std::cout << "inversedisplacementfieldfile: " << inversedisplacementfieldfile << std::endl;
    std::cout << "outputfile: \t\t\t" <<outputfile << std::endl;
  
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

    ImageToMeshFilterType::Pointer image2mesh = ImageToMeshFilterType::New();
    image2mesh->SetInput (reader->GetOutput());
    image2mesh->Update();

    WarperType::Pointer warper = WarperType::New();
    warper->SetDisplacementField (displacementfield);
    warper->SetInverseDisplacementField (inversedisplacementfield);
    warper->SetInput (image2mesh->GetOutput());
    warper->Update();

    TransformerType::Pointer transformer = TransformerType::New();
    transformer->SetTransform (transform);
    transformer->SetInput (warper->GetOutput());
    
    FilterType::Pointer gradientfilter = FilterType::New();
    gradientfilter->SetInput (0, transformer->GetOutput());
    gradientfilter->SetInput (1, transformer->GetOutput());
    gradientfilter->SetTransform (transform);
    gradientfilter->Update();
    

    MeshType::Pointer outputmesh1 = gradientfilter->GetOutput(0);
    MeshType::Pointer outputmesh2 = gradientfilter->GetOutput(1);
    MeshType::Pointer outputmesh3 = gradientfilter->GetOutput(2);

    
    TensorMeshIOType::Pointer meshwriter = TensorMeshIOType::New();
    meshwriter->SetInput (outputmesh1);
    meshwriter->SetFileName ("outputmesh.vtk");  
    meshwriter->Write();

    return EXIT_SUCCESS;
  
  }

}

  
