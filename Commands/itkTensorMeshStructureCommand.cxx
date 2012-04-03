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
#include <itkTensorMeshStructureCommand.h>

#include <itkTensorImageIO.h>
#include <itkTensorMeshIO.h>

#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "itkProlateSpheroidalStructureTensorMeshFilter.h"
#include "itkTensorImageToMeshFilter.h"
#include "itkTensorMeshToImageFilter.h"
#include "itkWarpTensorMeshFilter.h"
#include "itkProlateSpheroidalTransformTensorMeshFilter.h"
#include "itkLimitToAHAZoneImageFilter.h"
#include "GetPot.h"

namespace itk
{

  TensorMeshStructureCommand::TensorMeshStructureCommand()
  {
    m_ShortDescription = "Compute the structure tensor field of a tensor field in Prolate Spheroidal sense";
    m_LongDescription += m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-i    [input  tensor image]\n";
    m_LongDescription += "-pr   [prolate transform used]\n";
    m_LongDescription += "-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription += "-f2   [backward displacement field (default : backward.mha)]\n";  
    m_LongDescription += "-o    [output structure tensor field]\n";
    
  }

  TensorMeshStructureCommand::~TensorMeshStructureCommand()
  {}

  int TensorMeshStructureCommand::Execute (int narg, const char* arg[])
  {


    bool use_prolate = 1;
    
    
    typedef double                                      ScalarType;
    typedef TensorImageIO<ScalarType, 3, 3>             TensorIOType;
    typedef TensorMeshIO<ScalarType, 3, 3>              TensorMeshIOType;
    typedef TensorMeshIOType::TensorMeshType            MeshType;
    typedef Image<ScalarType, 3>                        ImageType;
    typedef ImageFileReader<ImageType>                  ImageReaderType;
    
    typedef ProlateSpheroidalStructureTensorMeshFilter<MeshType> FilterType;
    
    typedef FilterType::TensorType                      TensorType;
    typedef FilterType::TransformType                   TransformType;
    typedef TensorIOType::TensorImageType               TensorImageType;
    typedef itk::Vector<ScalarType, 3>                  DisplacementType;
    typedef itk::Image<DisplacementType, 3>             DisplacementFieldType;
    typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFileReaderType;
    typedef itk::TensorImageToMeshFilter<TensorType, 3> ImageToMeshFilterType;
    typedef itk::TensorMeshToImageFilter<TensorType, 3> MeshToImageFilterType;

    typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>  WarperType;
    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>   TransformerType;    

    typedef itk::LimitToAHAZoneImageFilter<TensorImageType>             AHALimiterType;    
    
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
    ImageReaderType::Pointer t_reader = ImageReaderType::New();
    
    reader->SetFileName(inputfile);
    t_reader->SetFileName(inputfile);
  
    try
    {
      reader->Read();
      t_reader->Update();
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

    TransformType::Pointer backtransform = TransformType::New();
    transform->GetInverse (backtransform);
    
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
    transformer->Update();
    
    MeshType::Pointer data;

    if (use_prolate)
      data = transformer->GetOutput();
    else
      data = image2mesh->GetOutput();
    
    data->DisconnectPipeline();

    
    FilterType::Pointer structurefilter = FilterType::New();
    structurefilter->SetInput (0, data);
    structurefilter->SetInput (1, data);
    structurefilter->SetTransform (transform);
    structurefilter->SetUsePiWorkAround(use_prolate);
    
    structurefilter->Update();

    MeshType::Pointer structure = structurefilter->GetOutput();
    structure->DisconnectPipeline();
    
    TransformerType::Pointer backtransformer = TransformerType::New();

    if (use_prolate)
    {
      backtransformer->SetTransform (backtransform);
      backtransformer->SetInput (structure);
      backtransformer->Update();
    }
    
    WarperType::Pointer backwarper = WarperType::New();
    if (use_prolate)
    {
      backwarper->SetDisplacementField (inversedisplacementfield);
      backwarper->SetInverseDisplacementField (displacementfield);
      backwarper->SetInput (backtransformer->GetOutput());
      backwarper->Update();    
    }

    MeshType::Pointer output;
    if (use_prolate)
      output = backwarper->GetOutput();
    else
      output = structure;
    
    output->DisconnectPipeline();

    MeshToImageFilterType::Pointer mesh2image = MeshToImageFilterType::New();
    mesh2image->SetInput (output);
    mesh2image->SetDomain (t_reader->GetOutput());
    mesh2image->Update();
    
    TensorIOType::Pointer outputwriter = TensorIOType::New();
    outputwriter->SetInput (mesh2image->GetOutput());
    outputwriter->SetFileName (outputfile);  
    outputwriter->Write();
    
    AHALimiterType::Pointer zonelimiter = AHALimiterType::New();
    zonelimiter->SetInput (mesh2image->GetOutput());
    zonelimiter->SetTransform (transform);
    zonelimiter->SetDisplacementField (displacementfield);
    zonelimiter->SetInverseDisplacementField (inversedisplacementfield);
    zonelimiter->CanineDivisionsOff();
    zonelimiter->SetAHASegmentationType (AHALimiterType::AHA_17_ZONES);
    
    MeshType::Pointer               zonestructure = MeshType::New();
    MeshType::PointsContainer::Pointer points     = MeshType::PointsContainer::New();
    MeshType::PointDataContainer::Pointer tensors = MeshType::PointDataContainer::New();
    zonestructure->SetPoints (points);
    zonestructure->SetPointData (tensors);
    points->Reserve (zonelimiter->GetNumberOfAHAZones());
    tensors->Reserve (zonelimiter->GetNumberOfAHAZones());
    
    for (unsigned int i=1; i<=zonelimiter->GetNumberOfAHAZones(); i++)
    {
      zonelimiter->SetAHAZone (i);
      ImageToMeshFilterType::Pointer zoneimage2mesh = ImageToMeshFilterType::New();
      zoneimage2mesh->SetInput (zonelimiter->GetOutput());
      zoneimage2mesh->Update();
      MeshType::Pointer zone = zoneimage2mesh->GetOutput();

      MeshType::PointType p = zonelimiter->GetZoneCentralPointCartesian();
      TensorType t (0.0);
      
      for (unsigned int j=0; j<zone->GetNumberOfPoints(); j++)
      {
	TensorType l (0.0); zone->GetPointData (j, &l);
	t += l.Log();
      }
      
      if (zone->GetNumberOfPoints())
      {
	t /= static_cast<ScalarType> (zone->GetNumberOfPoints());
	t = t.Exp();
      }
      
      zonestructure->SetPoint (i-1, p);
      zonestructure->SetPointData (i-1, t);
    }

    TensorMeshIOType::Pointer outputwriter2 = TensorMeshIOType::New();
    outputwriter2->SetInput (zonestructure);
    outputwriter2->SetFileName ("zones-structure.vtk");
    outputwriter2->Update();
    
    return EXIT_SUCCESS;
  
  }

}

  
