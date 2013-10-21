/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkReorderDWIsCommand.cxx 1 2010-05-21 14:00:33Z nt08 $
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
 * \file itkReorderDWIsCommand.cxx
 *
 * \author Nicolas Toussaint, King's College London, nicolas.toussaint@kcl.ac.uk
 * 
 * \url http://www.kcl.ac.uk/schools/medicine/research/imaging/
 *  
 */
#include <itkReorderDWIsCommand.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkNumericTraits.h>
#include <itksys/SystemTools.hxx>

#include <itkGradientFileReader.h>
#include <itkGradientFileWriter.h>

#include <sstream>
#include <fstream>
#include <iostream>
#include "GetPot.h"


namespace itk
{

  template<typename GradientListType>
  void EstimateOutputGradientList(GradientListType inputlist, bool replace, GradientListType &outputlist, std::vector<int> &ids)
  {
    int first_b0_id      = -1; // it is ALWAYS the case
    int last_b0_id       = -1; // it is ALWAYS the case
    int first_meandwi_id = -1; // it is NOT ALWAYS the case
    int last_meandwi_id  = -1; // it is NOT ALWAYS the case
  
    for (unsigned int i=0; i<inputlist.size(); i++)
    {
      if (!itk::NumericTraits<double>::IsPositive (inputlist[i].GetNorm()))
      {
	first_b0_id = i;
	break;
      }
    }
    if (first_b0_id == -1)
    {
      std::cerr<<"cannot find any b-0 in image!"<<std::endl;
      std::cerr<<"this is too weird, exiting"<<std::endl;
      std::exit (EXIT_FAILURE);
    }
    else
    {  
      for (int i=first_b0_id; i<(int)inputlist.size(); i++)
      {
	if (itk::NumericTraits<double>::IsPositive (inputlist[i].GetNorm()))
	{
	  last_b0_id = i-1;
	  break;
	}
      }
    }
  
    for (int i=last_b0_id+1; i<(int)inputlist.size(); i++)
    {
      if (!itk::NumericTraits<double>::IsPositive (inputlist[i].GetNorm()))
      {
	first_meandwi_id = i;
	break;
      }
    }
    if (first_meandwi_id == -1)
    {
      std::cerr<<"cannot find any mean diffusivity in image!"<<std::endl;
      std::cerr<<"we skip but do not exit"<<std::endl;
    }
    else
    {  
      for (int i=first_meandwi_id; i<(int)inputlist.size(); i++)
      {
	if (itk::NumericTraits<double>::IsPositive (inputlist[i].GetNorm()))
	{
	  last_meandwi_id = i-1;
	  break;
	}
      }
      if (last_meandwi_id == -1)
      {
	last_meandwi_id = inputlist.size()-1;
      }
    }

    
    std::cout<<"pattern found : "<<std::endl;
    std::cout<<"first b-0 at : "<<first_b0_id<<std::endl;
    std::cout<<"last b-0 at : "<<last_b0_id<<std::endl;
    std::cout<<"first mean diff. at : "<<first_meandwi_id<<std::endl;
    std::cout<<"last mean diff. at : "<<last_meandwi_id<<std::endl;


    outputlist.push_back (inputlist[first_b0_id]);
    ids.push_back (first_b0_id);
  
    std::cout<<"b-0 ids now inserted, size is : "<<ids.size()<<std::endl;
    for (unsigned int i=0; i<inputlist.size(); i++)
    {
      if (itk::NumericTraits<double>::IsPositive (inputlist[i].GetNorm()))
      {
	outputlist.push_back (inputlist[i]);
	ids.push_back (i);
      }
    }

    std::cout<<"DWIs ids now inserted, final size : "<<ids.size()<<std::endl;

  }

  template<typename Image4DType, typename ImageType, typename GradientListType>
  void ExtractB0Image(typename Image4DType::Pointer input, GradientListType gradientlist, typename ImageType::Pointer output)
  {
    typename Image4DType::SizeType           size = input->GetLargestPossibleRegion().GetSize();
    typename Image4DType::PointType        origin = input->GetOrigin();
    typename Image4DType::SpacingType     spacing = input->GetSpacing();
    typename Image4DType::DirectionType direction = input->GetDirection();

    typename ImageType::RegionType mregion;
    typename ImageType::SizeType msize;
    typename ImageType::PointType morigin;
    typename ImageType::SpacingType mspacing;
    typename ImageType::DirectionType mdirection;  
    for (unsigned int i=0; i<3; i++)
    {
      msize[i] = size[i];
      morigin[i] = origin[i];
      mspacing[i] = spacing[i];
      for (unsigned int j=0; j<3; j++)
	mdirection[i][j] = direction[i][j];
    }
    
    mregion.SetSize (msize);
    output->SetRegions (mregion);
    output->SetOrigin(morigin);
    output->SetSpacing(mspacing);  
    output->SetDirection(mdirection);
    output->Allocate();
    output->FillBuffer(static_cast<typename ImageType::PixelType>(0.0));
    
    typedef ExtractImageFilter< Image4DType, ImageType > ExtractorType;
    typename ExtractorType::Pointer  extractor  = ExtractorType::New();
    typename Image4DType::SizeType   buffersize = input->GetLargestPossibleRegion().GetSize(); buffersize[3] = 0;
    typename Image4DType::IndexType  bufferstart;  bufferstart.Fill (0);
    typename Image4DType::RegionType bufferregion; bufferregion.SetSize (buffersize); bufferregion.SetIndex (bufferstart);
    extractor->SetInput (input);
    extractor->SetDirectionCollapseToSubmatrix();
    
    unsigned counter = 0;
    for (unsigned int i = 0; i < gradientlist.size(); i++)
    {
      if (itk::NumericTraits<double>::IsPositive (gradientlist[i].GetNorm()))
	continue;
      bufferstart[3] = i;
      bufferregion.SetIndex (bufferstart);
      extractor->SetExtractionRegion (bufferregion);
      extractor->Update();
      typename ImageType::Pointer buffer = extractor->GetOutput();
      typename itk::ImageRegionIterator<ImageType> itIn(buffer, buffer->GetLargestPossibleRegion());
      typename itk::ImageRegionIterator<ImageType> itB0(output, output->GetLargestPossibleRegion());    
      while (!itB0.IsAtEnd())
      {
	itB0.Set (itB0.Get() + itIn.Get());
	++itIn; ++itB0;
      }
      counter++;
    }
    typename itk::ImageRegionIterator<ImageType> itB0(output, output->GetLargestPossibleRegion());    
    while (!itB0.IsAtEnd())
    {
      itB0.Set (itB0.Get() / (double)(counter));
      ++itB0;
    }
  }
  template <typename Image4DType, typename ImageType, typename GradientListType>
  void EstimateArithmeticcMean ( typename Image4DType::Pointer input, typename ImageType::Pointer b0, GradientListType gradientlist, double factor, typename ImageType::Pointer mean)
  {
    std::cout<<"arithmetic mean with a = "<<factor<<std::endl;
    
    mean->SetRegions (b0->GetLargestPossibleRegion());
    mean->SetOrigin(b0->GetOrigin());
    mean->SetSpacing(b0->GetSpacing());  
    mean->SetDirection(b0->GetDirection());
    mean->Allocate();
    mean->FillBuffer(static_cast<typename ImageType::PixelType>(0.0));

    typedef ExtractImageFilter< Image4DType, ImageType > ExtractorType;
    typename ExtractorType::Pointer  extractor  = ExtractorType::New();
    typename Image4DType::SizeType   buffersize = input->GetLargestPossibleRegion().GetSize(); buffersize[3] = 0;
    typename Image4DType::IndexType  bufferstart;  bufferstart.Fill (0);
    typename Image4DType::RegionType bufferregion; bufferregion.SetSize (buffersize); bufferregion.SetIndex (bufferstart);
    extractor->SetInput (input);
    extractor->SetDirectionCollapseToSubmatrix();
    
    unsigned int numberofdwis = 0;

    // 1) ACCUMULATE DWI SIGNAL FOR ARITHMETIC MEAN CALCULATION
    for (unsigned int i = 0; i < gradientlist.size(); i++)
    {
      if (!itk::NumericTraits<double>::IsPositive (gradientlist[i].GetNorm()))
	continue;

      numberofdwis++;
    
      bufferstart[3] = i;
      bufferregion.SetIndex (bufferstart);
      extractor->SetExtractionRegion (bufferregion);
      try
      {
	extractor->Update();
      }
      catch(itk::ExceptionObject &e)
      {
	std::cerr << e;
	exit(EXIT_FAILURE);
      }
      typename ImageType::Pointer buffer = extractor->GetOutput();
      
      typename itk::ImageRegionIterator<ImageType> itIn(buffer, buffer->GetLargestPossibleRegion());
      typename itk::ImageRegionIterator<ImageType> itMean(mean, mean->GetLargestPossibleRegion());
      
      while (!itIn.IsAtEnd())
      {
	// accumulating the pixel values (geometric mean)
	itMean.Set (itMean.Get() + itIn.Get());
	// iterate
	++itIn;
	++itMean;
      }
    }
    
    typename itk::ImageRegionIterator<ImageType> itMean(mean, mean->GetLargestPossibleRegion());

    // 2) ACTUAL DWI ARITHMETIC MEAN CALCULATION
    while (!itMean.IsAtEnd())
    {
      itMean.Set (factor * itMean.Get() / (double)(numberofdwis));      
      ++itMean;
    }
  }


  template <typename Image4DType, typename ImageType, typename GradientListType>
  void EstimateGeometricMean ( typename Image4DType::Pointer input, typename ImageType::Pointer b0, GradientListType gradientlist, double factor, typename ImageType::Pointer mean)
  {
    std::cout<<"geometric mean with f = "<<factor<<std::endl;

    mean->SetRegions (b0->GetLargestPossibleRegion());
    mean->SetOrigin(b0->GetOrigin());
    mean->SetSpacing(b0->GetSpacing());  
    mean->SetDirection(b0->GetDirection());
    mean->Allocate();
    mean->FillBuffer(static_cast<typename ImageType::PixelType>(1.0));
    
    typename ImageType::Pointer geometriccoefficientimage = ImageType::New();
    geometriccoefficientimage->SetRegions (b0->GetLargestPossibleRegion());
    geometriccoefficientimage->SetOrigin(b0->GetOrigin());
    geometriccoefficientimage->SetSpacing(b0->GetSpacing());  
    geometriccoefficientimage->SetDirection(b0->GetDirection());
    geometriccoefficientimage->Allocate();
    geometriccoefficientimage->FillBuffer(static_cast<typename ImageType::PixelType>(0.0));

    typedef ExtractImageFilter< Image4DType, ImageType > ExtractorType;
    typename ExtractorType::Pointer  extractor  = ExtractorType::New();
    typename Image4DType::SizeType   buffersize = input->GetLargestPossibleRegion().GetSize(); buffersize[3] = 0;
    typename Image4DType::IndexType  bufferstart;  bufferstart.Fill (0);
    typename Image4DType::RegionType bufferregion; bufferregion.SetSize (buffersize); bufferregion.SetIndex (bufferstart);
    extractor->SetInput (input);
    extractor->SetDirectionCollapseToSubmatrix();
    
    unsigned int numberofdwis = 0;

    // 1) ACCUMULATE DWI SIGNAL FOR GEOMETRIC MEAN CALCULATION AND GEOMETRIC COEFFICIENT CALCULATION
    for (unsigned int i = 0; i < gradientlist.size(); i++)
    {
      if (!itk::NumericTraits<double>::IsPositive (gradientlist[i].GetNorm()))
	continue;

      numberofdwis++;
    
      bufferstart[3] = i;
      bufferregion.SetIndex (bufferstart);
      extractor->SetExtractionRegion (bufferregion);
      try
      {
	extractor->Update();
      }
      catch(itk::ExceptionObject &e)
      {
	std::cerr << e;
	exit(EXIT_FAILURE);
      }
      typename ImageType::Pointer buffer = extractor->GetOutput();
      
      typename itk::ImageRegionIterator<ImageType> itIn(buffer, buffer->GetLargestPossibleRegion());
      typename itk::ImageRegionIterator<ImageType> itMean(mean, mean->GetLargestPossibleRegion());
      typename itk::ImageRegionIterator<ImageType> itGeom(geometriccoefficientimage, geometriccoefficientimage->GetLargestPossibleRegion());
      typename itk::ImageRegionIterator<ImageType> itB0(b0, b0->GetLargestPossibleRegion());
      
      while (!itIn.IsAtEnd())
      {
	// accumulating the pixel values (geometric mean)
	itMean.Set (itMean.Get() * itIn.Get());
	
	if ( (itB0.Get() > 0.0) && (itIn.Get() > 0.0) )
	  itGeom.Set (itGeom.Get() - std::log (itIn.Get() / itB0.Get()));
	// iterate
	++itIn;
	++itMean;
	++itGeom;
	++itB0;
      }
    }
    
    double expectedvalue = 0;
    unsigned int counter = 0;
    typename itk::ImageRegionIterator<ImageType> itMean(mean, mean->GetLargestPossibleRegion());
    typename itk::ImageRegionIterator<ImageType> itGeom(geometriccoefficientimage, geometriccoefficientimage->GetLargestPossibleRegion());
    typename itk::ImageRegionIterator<ImageType> itB0(b0, b0->GetLargestPossibleRegion());
    
    // 2) ACTUAL DWI GEOMETRIC MEAN CALCULATION AND GEOMETRIC COEFFICIENT EXPECTED VALUE
    while (!itMean.IsAtEnd())
    {
      typename ImageType::PixelType geometricmean = std::pow (itMean.Get(), 1.0 / (double)(numberofdwis));
      itMean.Set (geometricmean);
      itGeom.Set (itGeom.Get() / (double)(numberofdwis));
      
      if ( (itGeom.Get() > 0) && (itGeom.Get() < std::sqrt (2.0) ) )
      {
	counter++;
	expectedvalue += itGeom.Get();
      }
    
      // iterate
      ++itMean;
      ++itGeom;
    }
    
    expectedvalue /= (double)(counter);
    
    // 3) CORRECTION OF THE GEOMETRIC MEAN BY FACTOR AND EXPECTED VALUE
    itMean.GoToBegin();
    itGeom.GoToBegin();
    itB0.GoToBegin();
    
    while (!itMean.IsAtEnd())
    {
      itMean.Set ( itMean.Get() * std::exp (factor * expectedvalue + itGeom.Get()) );
	
      // iterate
      ++itMean;
      ++itGeom;
      ++itB0;
    }

    
    // END OF MEAN DIFFUSIVITY CALCULUS
    std::cout<<"calculating the actual mean diffusivity: done. mean geometric coefficient is "<<expectedvalue<<std::endl;  
    
    std::cout << "Writing coeff image to coeff.mha" << std::endl;
    typedef ImageFileWriter <ImageType> Image3DWriterType;  
    typename Image3DWriterType::Pointer coeffwriter = Image3DWriterType::New();
    coeffwriter->SetFileName ("coeff.mha");
    coeffwriter->SetInput(geometriccoefficientimage);
    coeffwriter->Update();
  }


  template <typename Image4DType, typename ImageType, typename GradientListType>
  void FillOutput( typename Image4DType::Pointer input, typename ImageType::Pointer firstimage, GradientListType gradientlist, std::vector<int> ids, typename Image4DType::Pointer output)
  {
    typename Image4DType::SizeType   size = input->GetLargestPossibleRegion().GetSize(); size[3] = ids.size();
    typename Image4DType::RegionType region; region.SetSize (size);
    
    output->SetRegions (region);
    output->SetOrigin(input->GetOrigin());
    output->SetSpacing(input->GetSpacing());  
    output->SetDirection(input->GetDirection());
    output->Allocate();
    output->FillBuffer(static_cast<typename ImageType::PixelType>(1.0));

    typedef ExtractImageFilter< Image4DType, ImageType > ExtractorType;
    typename ExtractorType::Pointer  extractor  = ExtractorType::New();
    typename Image4DType::SizeType   buffersize = input->GetLargestPossibleRegion().GetSize(); buffersize[3] = 0;
    typename Image4DType::IndexType  bufferstart;  bufferstart.Fill (0);
    typename Image4DType::RegionType bufferregion; bufferregion.SetSize (buffersize); bufferregion.SetIndex (bufferstart);
    extractor->SetInput (input);
    extractor->SetDirectionCollapseToSubmatrix();

    std::cout << "Filling output with extracted DWI..."<< std::endl;
    typename itk::ImageRegionIterator<Image4DType> itOut(output, output->GetLargestPossibleRegion());
  
    for (unsigned int i = 0; i < ids.size(); i++)
    {
      int id = ids[i];
  
      bufferstart[3] = id;
      bufferregion.SetIndex (bufferstart);
      extractor->SetExtractionRegion (bufferregion);
      try
      {
	extractor->Update();
      }
      catch(itk::ExceptionObject &e)
      {
	std::cerr << e;
	exit(EXIT_FAILURE);
      }

      typename ImageType::Pointer buffer = extractor->GetOutput();

      typename itk::ImageRegionIterator<ImageType> itIn(buffer, buffer->GetLargestPossibleRegion());
      typename itk::ImageRegionIterator<ImageType> itFirst(firstimage, firstimage->GetLargestPossibleRegion());
      // itGeom.GoToBegin();
    
      while (!itIn.IsAtEnd())
      {
	if (i == 0)
	  itOut.Set (itFirst.Get());
	else
	  itOut.Set (itIn.Get());

	++itIn;
	++itOut;
	++itFirst;
      }
    }

    std::cout << "Filling output with extracted DWI: Done." << std::endl;
  }
  
  
  ReorderDWIsCommand::ReorderDWIsCommand()
  {
    m_ShortDescription = "Reorders DWIs in a coherent manner starting with the b-0";
    m_LongDescription = m_ShortDescription;
    m_LongDescription = "\n\nUsage:\n";
    m_LongDescription +="-i [input image]\n";
    m_LongDescription +="-g [gradient file]\n";
    m_LongDescription +="-t [0/1/2] (replace B-0 with calculated mean-diffusivity )\n";
    m_LongDescription +="-a arithmetic factor (1.0 is unchanged)\n";
    m_LongDescription +="-f geometric factor (0.0 is unchanged)\n";
    m_LongDescription +="-o [output file base (.mha and .grad used])\n";
  }

  ReorderDWIsCommand::~ReorderDWIsCommand()
  {}

  int ReorderDWIsCommand::Execute (int narg, const char* arg[])
  {

    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return EXIT_FAILURE;
    }
    
    const bool IsInputPresent        = cl.search(2,"-i","-I");
    const bool IsOutputPresent       = cl.search(2,"-o","-O");
    const bool IsGradientsPresent    = cl.search(2,"-g","-G");
    if(!IsInputPresent || !IsGradientsPresent || !IsOutputPresent )
    {
      std::cerr << "Error: Input and (or) output not set." << std::endl;
      return EXIT_FAILURE;
    }

    const char*        fileIn          = cl.follow("NoFile",2,"-i","-I");
    const char*        fileGrad        = cl.follow("NoFile",2,"-g","-G");
    const char*        fileOut         = cl.follow("NoFile",2,"-o","-O");
    const unsigned int replacementtype = cl.follow(0, 2,"-t","-T");
    const double       afactor         = cl.follow(1.0, 2, "-a", "-A");
    const double       gfactor         = cl.follow(0.0, 2, "-f", "-F");

    typedef double                                            ScalarType;  
    typedef itk::Image<ScalarType, 4>                         Image4DType;
    typedef itk::Image<ScalarType, 3>                         ImageType;
    typedef itk::ImageFileReader <Image4DType>                ImageReaderType;  
    typedef itk::ImageFileWriter <Image4DType>                ImageWriterType;  
    typedef itk::ImageFileWriter <ImageType>                  Image3DWriterType;  
    typedef itk::ExtractImageFilter< Image4DType, ImageType > ImageExtractorType;
    
    ImageReaderType::Pointer    imagereader    = ImageReaderType::New();
    ImageWriterType::Pointer    imagewriter    = ImageWriterType::New();
    ImageExtractorType::Pointer imageextractor = ImageExtractorType::New();
    
    std::cout<<"reading input file"<<std::endl;  
    imagereader->SetFileName (fileIn);
    try
    {
      imagereader->Update();
    } catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return EXIT_FAILURE;
    }
    std::cout<<"done."<<std::endl;
    Image4DType::Pointer inputimage = imagereader->GetOutput();
  
    std::cout<<"reading grad file"<<std::endl;
    itk::GradientFileReader::Pointer gradientreader = itk::GradientFileReader::New();
    gradientreader->SetFileName (fileGrad);
    try
    {
      gradientreader->Update();
    } catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return EXIT_FAILURE;
    }  
    std::cout<<"done."<<std::endl;
    itk::GradientFileReader::VectorListType inputgradients = gradientreader->GetGradientList();  


    // INSTANTIATE
    itk::GradientFileReader::VectorListType outputgradients;
    std::vector<int> ids;
    ImageType::Pointer b0image             = ImageType::New();
    ImageType::Pointer arithmeticmeanimage = ImageType::New();    
    ImageType::Pointer geometricmeanimage  = ImageType::New();
    
    try
    {
      
      // ESTIMATE THE OUTPUT GRADIENT LIST
      EstimateOutputGradientList<itk::GradientFileReader::VectorListType> (inputgradients, replacementtype, outputgradients, ids);
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e;
      exit(EXIT_FAILURE);
    }

    try
    {
      // EXTRACT THE B0 IMAGE
      ExtractB0Image<Image4DType,ImageType,itk::GradientFileReader::VectorListType> (inputimage, inputgradients, b0image);  
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e;
      exit(EXIT_FAILURE);
    }

    ImageType::Pointer firstimage;
    switch (replacementtype)
    {
	case 1:
	  try
	  {
	    // CALCULATE THE ARITHMETIC MEAN OF THE DIFFUSION MIMAGES
	    EstimateArithmeticcMean<Image4DType,ImageType,itk::GradientFileReader::VectorListType>(inputimage, b0image, inputgradients, afactor, arithmeticmeanimage);
	  }
	  catch(itk::ExceptionObject &e)
	  {
	    std::cerr << e;
	    exit(EXIT_FAILURE);
	  }

    	  firstimage = arithmeticmeanimage;
	  break;
	case 2:
	  try
	  { 
	    // CALCULATE THE GEOMETRIC MEAN OF THE DIFFUSION MIMAGES
	    EstimateGeometricMean<Image4DType,ImageType,itk::GradientFileReader::VectorListType>(inputimage, b0image, inputgradients, gfactor, geometricmeanimage);
	  }
	  catch(itk::ExceptionObject &e)
	  {
	    std::cerr << e;
	    exit(EXIT_FAILURE);
	  }
	  
	  firstimage = geometricmeanimage;
	  break;
	case 0:
	default:
	  firstimage = b0image;
	  break;
    }    

    // FILL OUTPUT WITH THE RIGHT IMAGES
    Image4DType::Pointer output = Image4DType::New();
    FillOutput<Image4DType,ImageType,itk::GradientFileReader::VectorListType>(inputimage, firstimage, inputgradients, ids, output);

    // WRITE THE OUTPUTS
    std::ostringstream outputfilename;
    outputfilename << fileOut <<".mha";
    std::ostringstream outputgradientfilename;
    outputgradientfilename << fileOut <<".grad";
  
    std::cout << "Writing image to : "<< outputfilename.str().c_str () << std::endl;
    imagewriter->SetFileName (outputfilename.str().c_str ());
    imagewriter->SetInput(output);
    imagewriter->Update();

    std::cout << "Writing gradients to : "<< outputgradientfilename.str().c_str () << std::endl;
    itk::GradientFileWriter::Pointer gradientwriter = itk::GradientFileWriter::New();
    gradientwriter->SetFileName (outputgradientfilename.str().c_str ());
    gradientwriter->SetGradientList (outputgradients);
    gradientwriter->Update();
    
    return EXIT_SUCCESS;

  }
}
