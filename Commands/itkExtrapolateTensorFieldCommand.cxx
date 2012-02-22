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

  =========================================================================*/


// 1/ read sparse tensor field (e.g. vtkPointSet filled with tensor information)
// 1.bis/ inputs : desired spacing, lambda, sigma, stop-criterion, 3 points defining ellipsoid, etc
// 2/ Use the normalized Gaussian extrapolation scheme to have a dense tensor field to be used as initialization
// 3/ do a basis change of our tensor fieldSSSS to get them in prolate sph. systen (with the ellipsoid as input)
// 4/ use an extrapolation scheme with attachement to data (sparse) and regularization (laplacian) to get
// the final dense tensor field
// 5/ export the energy evolution

// after that : compare with the original dense data with a similarity measure
#include <itkExtrapolateTensorFieldCommand.h>

#include <itkImageFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include <itkBSplineKernelFunction.h>
#include <itkGaussianKernelFunction.h>
#include <itksys/SystemTools.hxx>

#include <itkTensorImageIO.h>
#include <itkTensorMeshIO.h>
#include <itkKaiserBesselKernelFunction.h>
#include <itkSparseTensorsDiffusionTensorImageFilter.h>
#include <itkExtrapolateTensorField.h>
#include <itkTensorMeshToImageFilter.h>

#include "GetPot.h"

namespace itk
{

  ExtrapolateTensorFieldCommand::ExtrapolateTensorFieldCommand()
  {
    m_ShortDescription = "\nTensor Field Dense Approximation from Sparse Data\n\n";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "Usage:\n";
    
    m_LongDescription +="-p     [use prolate spheroid (default : 0.0)]\n";
    m_LongDescription +="-i     [input tensor unstructured grid (default : input.vtk)]\n";
    m_LongDescription +="-d     [domain of diffusion (default : domain.mha)]\n";
    m_LongDescription +="-f1    [forward displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2    [backward displacement field (default : backward.mha)]\n";
    m_LongDescription +="-pr    [prolate used (default : prolate.lms)]\n";
    m_LongDescription +="-k     [interpolation kernel to use [0: Gaussian][1: B-Spline][2: Kaiser-Bessel](default: 0)]\n";
    m_LongDescription +="-o     [output tensor image file after gaussian initialization (default : gaussiantensors.mha)]\n";
    m_LongDescription +="-u     [kernel size list file to use]\n";
  }

  ExtrapolateTensorFieldCommand::~ExtrapolateTensorFieldCommand()
  {}


  int ExtrapolateTensorFieldCommand::Execute (int narg, const char* arg[])
  {
    
    bool force_reading_all = 1;
    
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    bool useprolatesystem    = cl.follow(false,2,"-p","-P");
    const char* inputfile      = NULL;
    
    inputfile = cl.follow("input.vtk",2,"-i","-I");
    
    const char* domainfile = cl.follow("domain.mha",2,"-d","-D");

    const char* prolatefile         = cl.follow("prolate.lms",2,"-pr","-PR");
    const char* displacementfieldfile = cl.follow("forward.mha",2,"-f1","-F1");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char* gaussiantensorfile = cl.follow("output.mha",2,"-o","-O");
    const unsigned int kerneltouse   = cl.follow(0, 2,"-k","-K");
    const char* kernelfile    = cl.follow("nofile", 2,"-u","-U");
  
    std::cout << "Processing tensor extrapolation with following arguments: " << std::endl;
    std::cout << "useprolatesystem: \t\t" << useprolatesystem << std::endl;
    std::cout << "prolatefile: \t\t\t" << prolatefile << std::endl;
    std::cout << "displacementfieldfile: \t\t" << displacementfieldfile << std::endl;
    std::cout << "inversedisplacementfieldfile: \t" << inversedisplacementfieldfile << std::endl;
    std::cout << "inputfile: \t\t\t" << inputfile << std::endl;
    std::cout << "domainfile: \t\t\t" << domainfile << std::endl;
    std::cout << "outputfile: \t\t\t" <<gaussiantensorfile << std::endl;
    std::cout << "kerneltouse: \t\t" <<kerneltouse << std::endl;
    std::cout << std::flush;

  
    // typedefs
    typedef double ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3> TensorImageIOType;
    typedef itk::TensorMeshIO<ScalarType, 3, 3> TensorMeshIOType;

    typedef itk::GaussianKernelFunction GaussianKernelFunctionType;
    typedef itk::BSplineKernelFunction<3> BSplineKernelFunctionType;
    typedef itk::KaiserBesselKernelFunction KaiserBesselKernelType;

    typedef itk::ExtrapolateTensorField<ScalarType, 3>        ExtrapolateTensorFieldType;
    typedef ExtrapolateTensorFieldType::TransformType         TransformType;
    typedef ExtrapolateTensorFieldType::MeshType              MeshType;
    typedef ExtrapolateTensorFieldType::ImageType             ImageType;
    typedef ExtrapolateTensorFieldType::DisplacementFieldType DisplacementFieldType;
    typedef itk::ImageFileReader<ImageType>                   ImageFileReaderType;
    typedef itk::ImageFileReader<DisplacementFieldType>       DisplacementFileReaderType;
    typedef itk::TensorMeshToImageFilter<ExtrapolateTensorFieldType::TensorType, 3> TensorMeshToImageFilterType;



    std::cout<<"reading kernel list : "<<kernelfile<<std::endl;
    std::ifstream inputliststream (kernelfile);
    if(inputliststream.fail())
    {
      std::cerr << "Unable to open file: " << kernelfile << std::endl;
      std::exit (EXIT_FAILURE);
    }
    unsigned int NumberOfKernels = 0;
    inputliststream >> NumberOfKernels;
    std::cout<<"encountered : "<<NumberOfKernels<<std::endl;
  
    std::string sline = "";
    itksys::SystemTools::GetLineFromStream(inputliststream, sline);
  
    std::vector<double*> kernellist;
    for (unsigned int N=0; N<NumberOfKernels; N++)
    {
      std::string line = "";
      itksys::SystemTools::GetLineFromStream(inputliststream, line);
      itksys_ios::istringstream parse ( line );
      double* kernel = new double[3];
      for (unsigned int i=0; i<3; i++)
	parse >> kernel[i];
    
      kernellist.push_back (kernel);
    }

    std::cout<<"found "<<NumberOfKernels<<" kernels :"<<std::endl;
    double* alphas = new double[3*kernellist.size()];
    for (unsigned int N=0; N<NumberOfKernels; N++)
    {
      std::cout<<N<<": "
	       <<kernellist[N][0]<<" : "
	       <<kernellist[N][1]<<" : "
	       <<kernellist[N][2]<<std::endl;
      alphas[3*N+0] = 1.0 * kernellist[N][0];
      alphas[3*N+1] = 1.0 * kernellist[N][1];
      alphas[3*N+2] = 1.0 * kernellist[N][2];
    
    }
  
    // instantiation
    // read the domain image
    std::cout << "Reading domain: " << domainfile << std::flush;
    ImageFileReaderType::Pointer domainreader = ImageFileReaderType::New();
    domainreader->SetFileName(domainfile);
    try
    {
      domainreader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << " Done." << std::endl;
    ImageType::Pointer domain = domainreader->GetOutput();
  
    std::cout << "Reading input tensors: " << inputfile << std::endl;
    TensorMeshIOType::Pointer reader = TensorMeshIOType::New();
    reader->SetFileName (inputfile);
    reader->Read();
    MeshType::Pointer data = reader->GetOutput();
    std::cout<<" Done."<<std::endl;

    TransformType::Pointer transform = NULL;
    DisplacementFieldType::Pointer displacementfield = NULL;
    DisplacementFieldType::Pointer inversedisplacementfield = NULL;
  
    if (useprolatesystem || force_reading_all)
    {
      std::cout<<"reading transform "<<prolatefile<<std::endl;
      itk::TransformFactory<TransformType>::RegisterTransform ();
      itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
      transformreader->SetFileName( prolatefile );
      transformreader->Update();
      transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
      TransformType::Pointer transform_inverse = TransformType::New();
      transform->GetInverse(transform_inverse);
      std::cout << " Done." << std::endl;
  
      std::cout << "Reading forward field: " << displacementfieldfile << std::flush;
      DisplacementFileReaderType::Pointer displacementreader1 = DisplacementFileReaderType::New();
      displacementreader1->SetFileName(displacementfieldfile);
      try
      {
	displacementreader1->Update();
      }
      catch(itk::ExceptionObject &e)
      {
	std::cerr << e << std::endl;
	return EXIT_FAILURE;
      }
      std::cout << " Done." << std::endl;
      displacementfield = displacementreader1->GetOutput();

      std::cout << "Reading backward field: " << inversedisplacementfieldfile << std::flush;
      DisplacementFileReaderType::Pointer displacementreader2 = DisplacementFileReaderType::New();
      displacementreader2->SetFileName(inversedisplacementfieldfile);
      try
      {
	displacementreader2->Update();
      }
      catch(itk::ExceptionObject &e)
      {
	std::cerr << e << std::endl;
	return EXIT_FAILURE;
      }
      std::cout << " Done." << std::endl;
      inversedisplacementfield = displacementreader2->GetOutput();
    }

    ExtrapolateTensorFieldType::Pointer extrapolator = ExtrapolateTensorFieldType::New();
    extrapolator->SetInput (data);
    extrapolator->SetDomain (domain);
    extrapolator->SetDisplacementField (displacementfield);
    extrapolator->SetInverseDisplacementField (inversedisplacementfield);
    extrapolator->SetTransform (transform);
    extrapolator->SetUseProlateCoordinates (useprolatesystem);
    extrapolator->SetAlpha (alphas);
  
    KaiserBesselKernelType::Pointer kaiserkernel = KaiserBesselKernelType::New();
    kaiserkernel->SetWindowSize (6.3);
    kaiserkernel->SetBeta (22.0);
  
    switch(kerneltouse)
    {
	case 1:
	  extrapolator->SetKernel (BSplineKernelFunctionType::New());
	  break;
	case 2:
	  extrapolator->SetKernel (kaiserkernel);
	  break;
	case 0:
	default:
	  extrapolator->SetKernel (GaussianKernelFunctionType::New());
	  break;
    }

    ExtrapolateTensorFieldType::InterpolatorType::LimiterType::Pointer limiter = extrapolator->GetInterpolator()->GetLimiter();
    limiter->CanineDivisionsOff();
    limiter->SetTransform (transform);
    limiter->SetDisplacementField (displacementfield);
    limiter->SetInverseDisplacementField (inversedisplacementfield);
    limiter->SetAHASegmentationType (ExtrapolateTensorFieldType::InterpolatorType::LimiterType::AHA_1_ZONE);
    limiter->CalculateZones();
  
    try
    {
      extrapolator->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << " done." << std::endl;
    ExtrapolateTensorFieldType::MeshType::Pointer output = extrapolator->GetOutput();
    output->DisconnectPipeline();
  
    std::cout << "Writing the output... "<<gaussiantensorfile<< " " << std::flush;
    std::string extension = itksys::SystemTools::GetFilenameLastExtension(gaussiantensorfile).c_str();
    if (strcmp (extension.c_str(), ".vtk") == 0)
    {
      TensorMeshIOType::Pointer writer = TensorMeshIOType::New();
      writer->SetFileName (gaussiantensorfile);
      writer->SetInput (output);
      writer->Write();
    }
    else
    {
      TensorMeshToImageFilterType::Pointer meshtoimage = TensorMeshToImageFilterType::New();
      meshtoimage->SetInput (output);
      meshtoimage->SetDomain (domain);
      meshtoimage->Update();
      TensorImageIOType::Pointer writer = TensorImageIOType::New();
      writer->SetFileName(gaussiantensorfile);
      writer->SetInput (meshtoimage->GetOutput());
      writer->Write();
    }
    std::cout<<" Done."<<std::endl;

    return EXIT_SUCCESS;

  }

  
}
