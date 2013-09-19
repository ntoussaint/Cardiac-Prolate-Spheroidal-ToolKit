/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkExtractKernelEnveloppe.cxx 1 2010-05-21 14:00:33Z nt08 $
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
 * \file itkExtractKernelEnveloppe.cxx
 * 
 * \author Nicolas Toussaint
 */
#include <itkExtractKernelsEnveloppeCommand.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itksys/SystemTools.hxx>

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include "GetPot.h"

#include <itkLimitToAHAZoneImageFilter.h>
#include <itkProlateDistanceToPointImageFilter.h>
#include <itkProlateSpheroidalTransform.h>
#include <itkResampleImageFilter.h>
#include <itkGaussianKernelFunction.h>

namespace itk
{

  ExtractKernelsEnveloppeCommand::ExtractKernelsEnveloppeCommand()
  {
    m_ShortDescription = "Extract distance image of kernels in Prolate Spheroidal coordinates";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-d    [input domain image (default : domain.mha)]\n";
    m_LongDescription +="-k    [kernel size file]\n";
    m_LongDescription +="-pr   [prolate transform used]\n";
    m_LongDescription +="-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription +="-o    [output image file]\n";
  }

  ExtractKernelsEnveloppeCommand::~ExtractKernelsEnveloppeCommand()
  {}

  int ExtractKernelsEnveloppeCommand::Execute (int narg, const char* arg[])
  {

    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout  << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char* domainfile                   = cl.follow("domain.mha",2,"-d","-D");
    const char* kernelfile                   = cl.follow("kernels.csv",2,"-k","-K");
    const char* prolatefile                  = cl.follow("prolate.pr",2,"-pr","-PR");
    const char* displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char* outputfile                   = cl.follow("output.vtk",2,"-o","-O");
    
    std::cout << "Processing kernel enveloppe extraction with following arguments: " << std::endl;
    std::cout << "domainfile: " << domainfile << std::endl;
    std::cout << "kernelfile: " << kernelfile << std::endl;
    std::cout << "prolatefile: " << prolatefile << std::endl;
    std::cout << "displacementfieldfile: " << displacementfieldfile << std::endl;
    std::cout << "inversedisplacementfieldfile: " << inversedisplacementfieldfile << std::endl;
    //std::cout << "zone: " << zone << std::endl;
    std::cout << "output: " << outputfile << std::endl;  
    std::cout << std::flush;
    
    // typedefs
    typedef double                                                         ScalarType;
    typedef itk::Image<ScalarType,3>                                       ImageType;
    typedef itk::ImageFileReader<ImageType>                                ImageFileReaderType;
    typedef itk::ImageFileWriter<ImageType>                                ImageFileWriterType;
    typedef itk::ImageRegionIterator<ImageType>                            ImageIteratorType;
    typedef itk::GradientImageFilter<ImageType>                            GradientImageFilterType;
    typedef GradientImageFilterType::OutputPixelType                       CovariantVectorType;
    typedef itk::Vector<ScalarType, 3>                                          VectorType;
    typedef itk::Image<VectorType, 3>                                      VectorImageType;
    typedef itk::Image<CovariantVectorType, 3>                             GradientImageType;
    typedef itk::Matrix<ScalarType, 3, 3>                                  MatrixType;
    typedef itk::LinearInterpolateImageFunction<ImageType, double>         InterpolatorType;
    typedef itk::Vector<ScalarType, 3>                                          DisplacementType;
    typedef itk::Image<DisplacementType, 3>                                DisplacementFieldType;
    typedef itk::ImageFileReader<DisplacementFieldType>                    DisplacementFileReaderType;

    typedef itk::ProlateSpheroidalTransform<ScalarType>                    TransformType;
    typedef itk::ProlateDistanceToPointImageFilter<ImageType>              DistanceMapperType;
    typedef itk::LimitToAHAZoneImageFilter<ImageType>                      AHALimiterType; 
    typedef TransformType::InputPointType                                  PointType;
    typedef itk::ResampleImageFilter<ImageType,ImageType>                  ResampleImageFilterType;

    typedef itk::GaussianKernelFunction<ScalarType> GaussianKernelFunctionType;

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
    for (unsigned int N=0; N<NumberOfKernels; N++)
    {
      std::cout<<N<<": "
	       <<kernellist[N][0]<<" : "
	       <<kernellist[N][1]<<" : "
	       <<kernellist[N][2]<<std::endl;
    }
  
    std::cout << "Reading input domain: " << domainfile <<"... "<<std::flush;
    ImageFileReaderType::Pointer imagereader = ImageFileReaderType::New();
    imagereader->SetFileName(domainfile);
    try
    {
      imagereader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ImageType::Pointer domain = imagereader->GetOutput();

    unsigned int resolutionfactor = 4;
    ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
    resampler->SetInput (imagereader->GetOutput());
    resampler->SetOutputOrigin (imagereader->GetOutput()->GetOrigin());
    resampler->SetOutputDirection (imagereader->GetOutput()->GetDirection());
    ImageType::SpacingType inputspacing = imagereader->GetOutput()->GetSpacing();
    ImageType::SpacingType outputspacing;
    ImageType::SizeType inputsize = imagereader->GetOutput()->GetLargestPossibleRegion().GetSize();
    ImageType::SizeType outputsize;
    for (unsigned int i=0; i<3; i++)
    {
      outputspacing[i] = inputspacing[i] / (double)(resolutionfactor);
      outputsize[i] = inputsize[i] * resolutionfactor;
    }
    resampler->SetOutputSpacing (outputspacing);
    resampler->SetSize (outputsize);
    try
    {
      resampler->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
  
    ImageType::Pointer domainfordistance = resampler->GetOutput();
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

  
    std::cout << "allocating..."<<std::endl;
  
    ImageType::DirectionType direction;
    direction.SetIdentity();
    ImageType::PointType origin;
    origin.Fill (0.0);
    ImageType::SpacingType spacing;
    spacing.Fill (1.0);
    ImageType::RegionType region;
    ImageType::SizeType size;
  
    for (unsigned int i=0; i<3; i++)
    {
      origin[i] = domainfordistance->GetOrigin()[i];
      spacing[i] = domainfordistance->GetSpacing()[i];
      size[i] = domainfordistance->GetLargestPossibleRegion().GetSize()[i];
    
      for (unsigned int j=0; j<3; j++)
	direction[i][j] = domainfordistance->GetDirection()[i][j];
    }
    region.SetSize (size);
  
    ImageType::Pointer outputimage = ImageType::New();
    outputimage->SetRegions (region);
    outputimage->SetOrigin(origin);
    outputimage->SetSpacing(spacing);  
    outputimage->SetDirection(direction);
    outputimage->Allocate();
    outputimage->FillBuffer (0.0);
  
    itk::ImageRegionIterator<ImageType> itOut(outputimage, outputimage->GetLargestPossibleRegion());

    double sqrt2pi = std::sqrt (2.0 * vnl_math::pi);
    GaussianKernelFunctionType::Pointer kernel = GaussianKernelFunctionType::New();
  
    for (unsigned int zone = 1; zone <= 17; zone++)
    {
      std::cout<<"finding the centralpointof zone: "<<zone<<"... "<<std::flush;
      AHALimiterType::Pointer zonelimiter = AHALimiterType::New();
      zonelimiter->SetInput (domain);
      zonelimiter->SetAHAZone (zone);
      zonelimiter->SetTransform (transform);
      zonelimiter->SetInverseDisplacementField (inversedisplacementfield);
      zonelimiter->SetDisplacementField (displacementfield);
      zonelimiter->Update();
      ImageType::PointType zonecenter = zonelimiter->GetZoneCentralPointCartesian();
      std::cout<<"Done."<<std::endl;

      std::cout<<"Estimating Prolate distance to zone "<<zone <<" central point..."<<std::flush;
      DistanceMapperType::Pointer DistanceMapper = DistanceMapperType::New();
      DistanceMapper->SetInput (domainfordistance);
      DistanceMapper->SetInverseDisplacementField (inversedisplacementfield);
      DistanceMapper->SetTransform (transform);
      DistanceMapper->SetPoint (zonecenter);

      double* alphas = kernellist[zone-1];
    
      DistanceMapper->SetAlpha (alphas);
    
      try
      {
	DistanceMapper->Update();
      }
      catch(itk::ExceptionObject &e)
      {
	std::cerr << e << std::endl;
	std::exit(EXIT_FAILURE);
      }
      std::cout << " Done." << std::endl;
    
      ImageType::Pointer distancemap = DistanceMapper->GetOutput();
    
      itk::ImageRegionIterator<ImageType> itIn(distancemap, distancemap->GetLargestPossibleRegion());
      itOut.GoToBegin();
    
      while (!itOut.IsAtEnd())
      {
	double formervalue = itOut.Get();
	double value = 0.0;
	if (formervalue > 0.000001)
	  value = std::max (sqrt2pi * kernel->Evaluate (itIn.Get()), formervalue);
	else
	  value = sqrt2pi * kernel->Evaluate (itIn.Get());
      
	itOut.Set (value);
	++itIn;
	++itOut;
      }
    }
  
    std::cout << "Writing distance map: " << outputfile <<"... "<<std::flush;
    ImageFileWriterType::Pointer imagewriter = ImageFileWriterType::New();
    imagewriter->SetFileName(outputfile);
    imagewriter->SetInput (outputimage);
    try
    {
      imagewriter->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << " Done." << std::endl;

    std::cout << "Writing bi-domain: bidomain.mha... "<<std::flush;
    ImageFileWriterType::Pointer imagewriter2 = ImageFileWriterType::New();
    imagewriter2->SetFileName("bidomain.mha");
    imagewriter2->SetInput (domainfordistance);
    try
    {
      imagewriter2->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << " Done." << std::endl;
  
  
    return EXIT_SUCCESS;
  }

}
