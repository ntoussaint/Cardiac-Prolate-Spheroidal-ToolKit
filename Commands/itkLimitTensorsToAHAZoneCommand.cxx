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
#include <itkLimitTensorsToAHAZoneCommand.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

#include <itkLimitToAHAZoneImageFilter.h>
#include <itkProlateSpheroidalTransformTensorMeshFilter.h>
#include <itkTensorImageIO.h>
#include <itkTensorMeshIO.h>
#include <itkWarpTensorMeshFilter.h>
#include <itkTensorImageToMeshFilter.h>

#include "GetPot.h"

namespace itk
{

  LimitTensorsToAHAZoneCommand::LimitTensorsToAHAZoneCommand()
  {
    m_ShortDescription = "Crop a a tensor image/mesh to an AHA zone in the Prolate Spheroidal sense";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input tensor mesh]\n";
    m_LongDescription +="-d    [input domain image]\n";
    m_LongDescription +="-f1   [displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2   [inverse displacement field (default : backward.mha)]\n";
    m_LongDescription +="-pr   [prolate transform used (default : prolate.pr)]\n";
    m_LongDescription +="-z    [AHA zone to take into account (1~17) (default: 0 : ALL ZONES)]\n";
    m_LongDescription +="-o    [output image zone file]\n";
  }

  LimitTensorsToAHAZoneCommand::~LimitTensorsToAHAZoneCommand()
  {}

  int LimitTensorsToAHAZoneCommand::Execute (int narg, const char* arg[])
  {
    
    typedef double                               ScalarType;
    typedef itk::TensorMeshIO<ScalarType, 3, 3>  TensorMeshIOType;
    typedef itk::TensorImageIO<ScalarType, 3, 3> TensorImageIOType;
    typedef TensorMeshIOType::TensorMeshType     MeshType;
    typedef MeshType::PixelType                  TensorType;
    typedef itk::Image<ScalarType,3>             ImageType;
    typedef itk::ImageFileReader<ImageType>      ImageReaderType;
    typedef itk::Vector<ScalarType, 3>           DisplacementType;
    typedef itk::Image<DisplacementType, 3>      DisplacementFieldType;
    typedef itk::ProlateSpheroidalTransform<ScalarType> TransformType;
    typedef itk::LimitToAHAZoneImageFilter<ImageType>   AHALimiterType;
    typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFileReaderType;
    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>   TransformerType;
    typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>  WarperType;
    typedef itk::TensorImageToMeshFilter<TensorType, 3>                 ImageToMeshFilterType;

    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    const char*  inputfile                    = cl.follow("input.vtk",2,"-i","-I");
    const char*  domainfile                   = cl.follow("domain.mha",2,"-d","-D");
    const char*  prolatefile                  = cl.follow("prolate.pr",2,"-pr","-PR");
    const char*  displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char*  inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char*  outputfile                   = cl.follow("output.vtk",2,"-o","-O");
    const unsigned int zone                   = cl.follow(0, 2,"-z","-Z");

    std::cout << "Processing AHA segments: " << std::endl;
    std::cout << "inputfile: \t\t\t" << inputfile << std::endl;
    std::cout << "domainfile: \t\t\t" << domainfile << std::endl;
    std::cout << "prolatefile: \t\t\t" << prolatefile << std::endl;
    std::cout << "displacementfieldfile: \t\t" << displacementfieldfile << std::endl;
    std::cout << "inversedisplacementfieldfile: \t\t" << inversedisplacementfieldfile << std::endl;
    std::cout << "outputfile: \t\t\t" <<outputfile << std::endl;
    std::cout << "zone:      \t\t\t" <<zone << std::endl;
  
    std::cout << std::flush;

    std::cout << "Reading input tensor image/mesh: " << inputfile <<"... "<< std::flush;
    std::string extension = itksys::SystemTools::GetFilenameLastExtension(inputfile).c_str();

    MeshType::Pointer input = NULL;
    if (strcmp (extension.c_str(), ".vtk") == 0)
    {
      TensorMeshIOType::Pointer reader = TensorMeshIOType::New();
      reader->SetFileName (inputfile);
      reader->Read();
      input = reader->GetOutput();
    }
    else
    {
      TensorImageIOType::Pointer reader = TensorImageIOType::New();
      reader->SetFileName(inputfile);    
      reader->Read();

      ImageToMeshFilterType::Pointer imagetomesh = ImageToMeshFilterType::New();
      imagetomesh->SetInput (reader->GetOutput());
      imagetomesh->Update();
      input = imagetomesh->GetOutput();
    }
    std::cout<<" Done."<<std::endl;
    
    std::cout << "Reading domain : " << domainfile <<"... "<< std::flush;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(domainfile);
    reader->Update();
    ImageType::Pointer domain = reader->GetOutput();
    std::cout << " Done." << std::endl;
    
    std::cout<<"reading prolate transform"<<std::endl;
    TransformType::Pointer transform = TransformType::New();  
    itk::TransformFactory<TransformType>::RegisterTransform ();
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
    transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    TransformType::Pointer transform_inverse = TransformType::New();
    transform->GetInverse(transform_inverse);
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

    std::cout << "Reading backward field: " << inversedisplacementfieldfile << std::flush;  
    DisplacementFileReaderType::Pointer    inversedisplacementreader = DisplacementFileReaderType::New();
    inversedisplacementreader->SetFileName(inversedisplacementfieldfile);
    try
    {
      inversedisplacementreader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << " Done." << std::endl;
    DisplacementFieldType::Pointer inversedisplacementfield = inversedisplacementreader->GetOutput();    
    
    std::cout<<"Warping..."<<std::endl;
    WarperType::Pointer warper1 = WarperType::New();
    warper1->SetInput (input);
    warper1->SetDisplacementField (displacementfield);
    warper1->SetInverseDisplacementField (inversedisplacementfield);
    warper1->Update();
    
    std::cout<<"Transforming..."<<std::endl;
    TransformerType::Pointer transformer1 = TransformerType::New();
    transformer1->SetInput (warper1->GetOutput());
    transformer1->SetTransform (transform);
    transformer1->Update();
    
    MeshType::Pointer data1 = transformer1->GetOutput();

    data1->DisconnectPipeline();

    std::cout<<"limiting the tensors "<<inputfile<<" to AHA zone: "<<zone<<std::endl;
    AHALimiterType::Pointer zonelimiter = AHALimiterType::New();
    zonelimiter->SetInput (domain);
    zonelimiter->SetAHAZone (zone);
    zonelimiter->SetTransform (transform);
    zonelimiter->SetDisplacementField (inversedisplacementfield);
    zonelimiter->SetInverseDisplacementField (inversedisplacementfield);
    zonelimiter->CanineDivisionsOff();
    zonelimiter->SetAHASegmentationType (AHALimiterType::AHA_17_ZONES);
    zonelimiter->SetVentricleSizes (15, 95 * vnl_math::pi / 180.0);
    zonelimiter->CalculateZones();

    MeshType::Pointer data2 = MeshType::New();
    MeshType::PointsContainer::Pointer pts    = MeshType::PointsContainer::New();
    MeshType::PointDataContainer::Pointer tns = MeshType::PointDataContainer::New();
    data2->SetPoints (pts);
    data2->SetPointData (tns);

    unsigned int counter = 0;
    
    for (unsigned int i=0; i<data1->GetNumberOfPoints(); i++)
    {
      MeshType::PointType p;
      TensorType s (0.0);
      p[0] = p[1] = p[2] = 0.0;
      data1->GetPoint     (i, &p);
      data1->GetPointData (i, &s);
      
      unsigned int z = zonelimiter->InWhichZoneIsPoint (p);
      if (!z)
      {
	std::cerr<<"point "<<p<<" does not belong to any zone : "<<z<<std::endl;
	continue;
      }
      if (z == zone)
      {
	data2->GetPoints()->InsertElement    (counter, p);
	data2->GetPointData()->InsertElement (counter, s);
	counter++;
      }
    }
    std::cout<<"Done."<<std::endl;
    
    std::cout<<"Transforming..."<<std::endl;
    TransformerType::Pointer transformer2 = TransformerType::New();
    transformer2->SetInput (data2);
    transformer2->SetTransform (transform_inverse);
    transformer2->Update();

    std::cout<<"Warping..."<<std::endl;
    WarperType::Pointer warper2 = WarperType::New();
    warper2->SetInput (transformer2->GetOutput());
    warper2->SetDisplacementField (inversedisplacementfield);
    warper2->SetInverseDisplacementField (displacementfield);
    warper2->Update();

    MeshType::Pointer output = warper2->GetOutput();
    output->DisconnectPipeline();
    
    std::cout<<"writing resuting zone to : "<<outputfile<<std::endl;
    TensorMeshIOType::Pointer writer = TensorMeshIOType::New();
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
