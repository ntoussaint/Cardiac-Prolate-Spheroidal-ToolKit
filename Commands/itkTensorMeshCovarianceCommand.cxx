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
#include <itkTensorMeshCovarianceCommand.h>

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

  TensorMeshCovarianceCommand::TensorMeshCovarianceCommand()
  {
    m_ShortDescription = "Compute the covariance tensor field of a tensor field in Prolate Spheroidal sense";
    m_LongDescription += m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-i    [input structure tensor image]\n";
    m_LongDescription += "-pr   [prolate transform used]\n";
    m_LongDescription += "-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription += "-f2   [backward displacement field (default : backward.mha)]\n";  
    m_LongDescription += "-o    [output structure tensor field]\n";
    
  }

  TensorMeshCovarianceCommand::~TensorMeshCovarianceCommand()
  {}

  int TensorMeshCovarianceCommand::Execute (int narg, const char* arg[])
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
    const char* inputfile                    = cl.follow("input.mha",2,"-i","-I");
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
  
    std::cout << "Reading input structure tensor image: " << inputfile <<"... "<< std::flush;
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

    std::cout<<"now try the covariance..."<<std::endl;
    
    std::cout<<"defining the limiter..."<<std::endl;
    AHALimiterType::Pointer zonelimiter = AHALimiterType::New();
    zonelimiter->SetVentricleSizes(15, 105 * vnl_math::pi / 180.0);
    
    zonelimiter->SetInput (reader->GetOutput());
    zonelimiter->SetTransform (transform);
    zonelimiter->SetDisplacementField (displacementfield);
    zonelimiter->SetInverseDisplacementField (inversedisplacementfield);
    zonelimiter->CanineDivisionsOff();
    zonelimiter->SetAHASegmentationType (AHALimiterType::AHA_17_ZONES);
    zonelimiter->CalculateZones();

    std::cout<<"allocating the final structure tensors..."<<std::endl;
    MeshType::Pointer               zonecovariance = MeshType::New();
    MeshType::PointsContainer::Pointer points      = MeshType::PointsContainer::New();
    MeshType::PointDataContainer::Pointer tensors  = MeshType::PointDataContainer::New();
    zonecovariance->SetPoints (points);
    zonecovariance->SetPointData (tensors);
    points->Reserve (zonelimiter->GetNumberOfAHAZones());
    tensors->Reserve (zonelimiter->GetNumberOfAHAZones());
    
    std::cout<<"allocating vector of zone structure (tensors) in "<<zonelimiter->GetNumberOfAHAZones()<<" zones..."<<std::endl;
    std::vector<MeshType::Pointer> zonestructures;
    for (unsigned int i=1; i<=zonelimiter->GetNumberOfAHAZones(); i++)
    {
      MeshType::Pointer s = MeshType::New();
      MeshType::PointsContainer::Pointer pts    = MeshType::PointsContainer::New();
      MeshType::PointDataContainer::Pointer tns = MeshType::PointDataContainer::New();
      s->SetPoints (pts);
      s->SetPointData (tns);
      zonestructures.push_back (s);
    }

    std::cout<<"filling zone structure (tensors)..."<<std::endl;
    for (unsigned int i=0; i<data->GetNumberOfPoints(); i++)
    {
      MeshType::PointType p;
      TensorType s (0.0);
      p[0] = p[1] = p[2] = 0.0;
      data->GetPoint (i, &p);
      data->GetPointData (i, &s);
      
      unsigned int z = zonelimiter->InWhichZoneIsPoint (p);
      if (!z)
      {
	std::cerr<<"point "<<p<<" does not belong to any zone : "<<z<<std::endl;
	continue;
      }
      
      unsigned int size = zonestructures[z-1]->GetNumberOfPoints();
      
      zonestructures[z-1]->GetPoints()->InsertElement (size, p);
      zonestructures[z-1]->GetPointData()->InsertElement (size, s);
    }

    std::cout<<"evaluating the covariance for each zone..."<<std::endl;
    for (unsigned int z=1; z<=zonestructures.size(); z++)
    {
      
      std::cout<<"====== zone "<<z<<" ====="<<std::endl;

      MeshType::Pointer zone = zonestructures[z-1];
      
      std::cout<<"number of structures in zone : "<<zone->GetNumberOfPoints()<<std::endl;
      
      MeshType::PointType meanpoint; meanpoint[0] = meanpoint[1] = meanpoint[2] = 0.0;
      DisplacementType meangradient; meangradient[0] = meangradient[1] = meangradient[2] = 0.0;
      
      std::vector<DisplacementType> gradients;
 
      for (unsigned int i=0; i<zone->GetNumberOfPoints(); i++)
      {
	MeshType::PointType p;
	zone->GetPoint (i,&p);
	TensorType s (0.0);
	zone->GetPointData (i,&s);

	DisplacementType gradient = s.GetEigenvector (2);

	gradients.push_back (gradient);

	meangradient   += gradient;
	for (unsigned int u=0; u<3; u++)
	  meanpoint[u] += p[u];
      }

      if (gradients.size())
	meangradient /= (double)(gradients.size());
      if (gradients.size())
	for (unsigned int u = 0; u<3; u++)
	  meanpoint[u] /= (double)(gradients.size());

      std::cout<<"meangradient : "<<meangradient<<std::endl;
      std::cout<<"meanpoint : "<<meanpoint<<std::endl;
      
      vnl_matrix_fixed<ScalarType,3,3> covariancematrix (0.0);
      
      for (unsigned int j=0; j<gradients.size(); j++)
      {
	DisplacementType gradient = gradients[j];
	
	for (unsigned int k=0; k<3; k++)
	  for (unsigned int l=0; l<3; l++)
	    covariancematrix[k][l] += (gradient[k])*(gradient[l]);
      }
      
      if (gradients.size())
	covariancematrix /= (double)(gradients.size());
      
      TensorType covariance (0.0);
      covariance.SetVnlMatrix (covariancematrix);

      covariance = covariance.Sqrt();
      covariance *= 10;
      
      zonecovariance->SetPoint (z-1, meanpoint);
      zonecovariance->SetPointData (z-1, covariance);
      
    }

    TransformType::Pointer inversetransform = TransformType::New();
    transform->GetInverse (inversetransform);

    TransformerType::Pointer transformer2 = TransformerType::New();
    transformer2->SetTransform (inversetransform);
    transformer2->SetInput (zonecovariance);
    transformer2->Update();

    WarperType::Pointer warper2 = WarperType::New();
    warper2->SetDisplacementField (inversedisplacementfield);
    warper2->SetInverseDisplacementField (displacementfield);
    warper2->SetInput (transformer2->GetOutput());
    warper2->Update();

    MeshType::Pointer outputmesh = warper2->GetOutput();
    outputmesh->DisconnectPipeline();
    
    std::cout<<"correcting the zone-centre position..."<<std::endl;
    for (unsigned int i=0; i < zonecovariance->GetNumberOfPoints(); i++)
    {
      zonelimiter->SetAHAZone (i+1);
      MeshType::PointType p = zonelimiter->GetZoneCentralPointCartesian();
      outputmesh->SetPoint (i, p);
    }
    
    std::cout<<"Writing zone covariance tensors to "<<outputfile<<"... "<<std::flush;
    TensorMeshIOType::Pointer outputwriter2 = TensorMeshIOType::New();
    outputwriter2->SetInput (outputmesh);
    outputwriter2->SetFileName (outputfile);
    outputwriter2->Write();
    std::cout<<"done."<<std::endl;
    
    return EXIT_SUCCESS;
  
  }

}

  
