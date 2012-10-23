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
    
    m_LongDescription += "-i    [input  tensor image/mesh]\n";
    m_LongDescription += "-d    [domain file]\n";
    m_LongDescription += "-pr   [prolate transform used]\n";
    m_LongDescription += "-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription += "-f2   [backward displacement field (default : backward.mha)]\n";  
    m_LongDescription += "-o    [output structure tensor field]\n";
    m_LongDescription += "-f    [gradient evaluation factor (default=1.0)]\n";
    m_LongDescription += "-th   [wall thickness in mm (default=12.0)]\n";
    m_LongDescription += "-ma   [wall max angle in deg. (default=95.0)]\n";
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
    const char* domainfile                   = cl.follow("nofile",2,"-d","-D");
    const char* displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char* outputfile                   = cl.follow("output.mha",2,"-o","-O");
    const double factor                      = cl.follow(1.0,2,"-f","-F");
    const double thickness                   = cl.follow(12,2,"-th","-th");
    const double maxangle                    = cl.follow(95,2,"-ma","-MA");
  
    std::cout << "inputfile: \t\t\t" << inputfile << std::endl;
    std::cout << "factor: \t\t\t" << factor << std::endl;
    std::cout << "prolatefile: " << prolatefile << std::endl;
    std::cout << "domainfile: " << domainfile << std::endl;
    std::cout << "displacementfieldfile: " << displacementfieldfile << std::endl;
    std::cout << "inversedisplacementfieldfile: " << inversedisplacementfieldfile << std::endl;
    std::cout << "outputfile: \t\t\t" <<outputfile << std::endl;
  
    std::cout << std::flush;
  
    std::cout << "Reading input tensor image/mesh: " << inputfile <<"... "<< std::flush;
    std::string extension = itksys::SystemTools::GetFilenameLastExtension(inputfile).c_str();

    MeshType::Pointer data = NULL;
    if (strcmp (extension.c_str(), ".vtk") == 0)
    {
      TensorMeshIOType::Pointer reader = TensorMeshIOType::New();
      reader->SetFileName (inputfile);
      reader->Read();
      data = reader->GetOutput();
    }
    else
    {
      TensorIOType::Pointer reader = TensorIOType::New();
      reader->SetFileName(inputfile);    
      reader->Read();

      ImageToMeshFilterType::Pointer imagetomesh = ImageToMeshFilterType::New();
      imagetomesh->SetInput (reader->GetOutput());
      imagetomesh->Update();
      data = imagetomesh->GetOutput();
    }
    std::cout<<" Done."<<std::endl;

    std::cout << "Reading domain : " << domainfile <<"... "<< std::flush;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(domainfile);
    reader->Update();
    ImageType::Pointer domain = reader->GetOutput();
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

    data->DisconnectPipeline();
    
    std::cout<<"switching coordinates... "<<std::flush;
    WarperType::Pointer warper = WarperType::New();
    warper->SetDisplacementField (displacementfield);
    warper->SetInverseDisplacementField (inversedisplacementfield);
    warper->SetInput (data);
    warper->Update();

    TransformerType::Pointer transformer = TransformerType::New();
    transformer->SetTransform (transform);
    transformer->SetInput (warper->GetOutput());
    transformer->Update();
    std::cout<<"done."<<std::endl;
    
    MeshType::Pointer dataforstructure;

    if (use_prolate)
      dataforstructure = transformer->GetOutput();
    else
      dataforstructure = data;
    
    dataforstructure->DisconnectPipeline();

    double position_max[3] = {-100,-100,-100};
    double position_min[3] = {+100,+100,+100};

    for (unsigned int i=0; i < dataforstructure->GetNumberOfPoints(); i++)
    {
      MeshType::PointType p;
      dataforstructure->GetPoint (i, &p);
      double h[3];
      transform->EvaluateScaleFactors (p.GetDataPointer(), h);
      
      for (unsigned int j=0; j<3; j++)
      {
	position_max[j] = std::max (h[j], position_max[j]);
	position_min[j] = std::min (h[j], position_min[j]);
      }
    }

    std::cout<<"bounds : "<<std::endl;
    std::cout << position_min[0] << " : " << position_max[0] <<std::endl;
    std::cout << position_min[1] << " : " << position_max[1] <<std::endl;
    std::cout << position_min[2] << " : " << position_max[2] <<std::endl;
    
    
    std::cout<<"evaluating structure... "<<std::flush;
    FilterType::Pointer structurefilter = FilterType::New();
    structurefilter->SetInput (0, dataforstructure);
    structurefilter->SetInput (1, dataforstructure);
    structurefilter->SetTransform (transform);
    structurefilter->SetUsePiWorkAround(use_prolate);
    structurefilter->SetGradientUFactor (factor);
    
    structurefilter->Update();
    std::cout<<"done."<<std::endl;

    MeshType::Pointer structure = structurefilter->GetOutput();
    structure->DisconnectPipeline();

    std::cout<<"switching coordinates... "<<std::flush;
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
    std::cout<<"done."<<std::endl;

    std::cout<<"transforming to image for AHA separation... "<<std::flush;
    MeshToImageFilterType::Pointer mesh2image = MeshToImageFilterType::New();
    mesh2image->SetInput (output);
    mesh2image->SetDomain (domain);
    mesh2image->Update();

    std::string outextension = itksys::SystemTools::GetFilenameLastExtension(outputfile).c_str();

    if (strcmp (outextension.c_str(), ".vtk") == 0)
    {
      TensorMeshIOType::Pointer writer = TensorMeshIOType::New();
      writer->SetInput (output);
      writer->SetFileName (outputfile);
      writer->Write();
    }
    else
    {  
      TensorIOType::Pointer outputwriter = TensorIOType::New();
      outputwriter->SetInput (mesh2image->GetOutput());
      outputwriter->SetFileName (outputfile);
      outputwriter->Write();
    }
    
    AHALimiterType::Pointer zonelimiter = AHALimiterType::New();
    zonelimiter->SetInput (mesh2image->GetOutput());
    zonelimiter->SetTransform (transform);
    zonelimiter->SetDisplacementField (displacementfield);
    zonelimiter->SetInverseDisplacementField (inversedisplacementfield);
    zonelimiter->CanineDivisionsOff();
    zonelimiter->SetVentricleSizes (thickness, maxangle * vnl_math::pi / 180.0);
    
    zonelimiter->SetAHASegmentationType (AHALimiterType::AHA_17_ZONES);
    std::cout<<"done."<<std::endl;

    std::cout<<"computing structure covariance within zones... "<<std::flush;
    MeshType::Pointer               zonestructure = MeshType::New();
    MeshType::PointsContainer::Pointer points     = MeshType::PointsContainer::New();
    MeshType::PointDataContainer::Pointer tensors = MeshType::PointDataContainer::New();
    zonestructure->SetPoints (points);
    zonestructure->SetPointData (tensors);
    points->Reserve (zonelimiter->GetNumberOfAHAZones());
    tensors->Reserve (zonelimiter->GetNumberOfAHAZones());
    
    for (unsigned int i=1; i<=zonelimiter->GetNumberOfAHAZones(); i++)
    {
      std::cout<<"======= zone "<<i<<" ======= "<<std::flush;
      zonelimiter->SetAHAZone (i);
      ImageToMeshFilterType::Pointer zoneimage2mesh = ImageToMeshFilterType::New();
      zoneimage2mesh->SetInput (zonelimiter->GetOutput());
      zoneimage2mesh->Update();
      MeshType::Pointer zone = zoneimage2mesh->GetOutput();

      MeshType::PointType p = zonelimiter->GetZoneCentralPointCartesian();
      TensorType t (0.0);
      
      double ratio = 1.3;

      std::cout<<"( N = "<<zone->GetNumberOfPoints()<<") "<<std::flush;
      for (unsigned int j=0; j<zone->GetNumberOfPoints(); j++)
      {
	TensorType l (0.0); zone->GetPointData (j, &l);
	TensorType::MatrixType U;
	TensorType::MatrixType D;
	D.Fill (0.0);
	D[0][0] = l.GetEigenvalue (0);
	D[1][1] = l.GetEigenvalue (1);
	D[2][2] = l.GetEigenvalue (2) * (ratio);
	for (unsigned int a=0; a<3; a++)
	  for (unsigned int b=0; b<3; b++)
	  U[a][b] = l.GetEigenvector (b)[a];
	l.SetVnlMatrix (U * D * U.GetTranspose());
      
	t += l.Log();
      }
      
      if (zone->GetNumberOfPoints())
      {
	t /= static_cast<ScalarType> (zone->GetNumberOfPoints());
	t = t.Exp();
      }

      zonestructure->SetPoint (i-1, p);
      zonestructure->SetPointData (i-1, t);
      std::cout<<"done."<<std::endl;
    }

    std::cout<<"Writing zone structure tensors to "<<"zones-structure.vtk"<<"... "<<std::flush;
    TensorMeshIOType::Pointer outputwriter2 = TensorMeshIOType::New();
    outputwriter2->SetInput (zonestructure);
    outputwriter2->SetFileName ("zones-structure.vtk");
    outputwriter2->Write();
    std::cout<<"done."<<std::endl;
    
    return EXIT_SUCCESS;
  
  }

}

  
