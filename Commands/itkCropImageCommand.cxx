/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkCropImageCommand.cxx 1 2010-05-21 14:00:33Z nt08 $
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
 * \file itkCropImageCommand.cxx
 *
 * \author Nicolas Toussaint, King's College London, nicolas.toussaint@kcl.ac.uk
 * 
 * \url http://www.kcl.ac.uk/schools/medicine/research/imaging/
 *  
 */
#include <itkCropImageCommand.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include <itkRegionOfInterestImageFilter.h>

#include <itkProlateSpheroidalTransform.h>

#include <itksys/SystemTools.hxx>
#include "GetPot.h"

namespace itk
{

  /**
     Finds the region around the centre of the ventricle
     with a size of diameter x diameter

     The centre of the ventricle is found as the intersection between
     the image plane passing through the origin, and the line passing
     through the ellipsoid focii.
     
   */
  template <typename ImageType, typename TransformType>
  typename ImageType::RegionType
  FindRegion ( typename ImageType::Pointer image, typename TransformType::Pointer transform, double diameter)
  {
    typedef typename TransformType::ScalarType ScalarType;
    typedef vnl_matrix< ScalarType >           InternalMatrixType;
    typedef vnl_svd< ScalarType >              SolverType;
    
    // Recovering information from the transform (the focii positions)    
    typename TransformType::PointType f1 = transform->GetFocus1();
    typename TransformType::PointType f2 = transform->GetFocus2();

    // Recovering information from the image
    typename TransformType::VectorType spacing = image->GetSpacing().GetVnlVector().extract (3).data_block();
    typename TransformType::PointType P = image->GetOrigin().GetVnlVector().extract (3).data_block();
    InternalMatrixType direction = image->GetDirection().GetVnlMatrix().extract (3,3);
    typename TransformType::VectorType U = direction.get_column(0).data_block();
    typename TransformType::VectorType V = direction.get_column (1).data_block();

    // definition of vectors
    typename TransformType::VectorType axis = f2 - f1;
    typename TransformType::VectorType PF1 = f1 - P;

    
    // Finding the intersection between the image plane and the focii line
    // can be done by parametrizing both conditions:
    // (1) centre = P + alpha[0] * U + alpha[1] * V;
    // (2) center = f1 - alpha[2] * axis;
    // therefore
    // (3) alpha[0] * U + alpha[1] * V + alpha[2] * axis = f1 - P;

    // We represent (3) as a matrix linear problem :
    
    InternalMatrixType A (3,3);
    InternalMatrixType B (3,1);
    A.set_identity(); B.set_identity();
    
    A.set_column (0, U.GetDataPointer());
    A.set_column (1, V.GetDataPointer());
    A.set_column (2, axis.GetDataPointer());
    B.set_column (0, PF1.GetDataPointer());
    
    /// solve A . alpha = B
    SolverType solver (A);
    InternalMatrixType alpha = solver.solve (B);
    
    typename TransformType::PointType centre = P + alpha.get (0,0) * U + alpha.get (1,0) * V;
        
    typename TransformType::PointType origin = centre - 0.5 * diameter * (U + V);
    typename ImageType::PointType origin4d;
    for (unsigned int i=0; i<3; i++) origin4d[i] = origin[i];
    // origin4d[3] = 0;
    typename ImageType::IndexType index;
    image->TransformPhysicalPointToIndex (origin4d, index);
    
    typename ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    for (unsigned int i=0; i<2; i++) size[i] = vnl_math_rnd (diameter / spacing[i]);
    
    typename ImageType::RegionType region;
    region.SetIndex (index);
    region.SetSize (size);

    return region;
  }
  
  CropImageCommand::CropImageCommand()
  {
    m_ShortDescription = "crop a 4D to a certain diameter around the central long axis";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input 4D DWI image (default : input.mha)]\n";    
    m_LongDescription +="-pr   [prolate transform]\n";
    m_LongDescription +="-d    [diameter to use (default is 85mm)]\n";
    m_LongDescription +="-o    [output cropped image]\n";
  }

  CropImageCommand::~CropImageCommand()
  {}

  int CropImageCommand::Execute (int narg, const char* arg[])
  {
    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout  << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char* inputfile   = cl.follow("input.mha",2,"-i","-I");
    const double diameter   = cl.follow(85, 2,"-d","-D");
    const char* prolatefile = cl.follow("prolate.tr",2,"-pr","-PR");
    const char* outputfile  = cl.follow("output.csv",2,"-o","-O");

    typedef double                             ScalarType;
    typedef ProlateSpheroidalTransform<double> TransformType;
    
    std::cout<<"reading transform "<<prolatefile<<"..."<<std::endl;
    TransformFileReader::Pointer transformreader = TransformFileReader::New();
    itk::TransformFactory<TransformType>::RegisterTransform ();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
    std::cout << " Done." << std::endl;
    TransformType::Pointer transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    
    // first read the image component type
    typedef itk::Image<short,4> d_ImageType;
    typedef itk::ImageFileReader<d_ImageType> d_ReaderType;
    d_ReaderType::Pointer d_reader = d_ReaderType::New();
    d_reader->SetFileName(inputfile);
    d_reader->GenerateOutputInformation();

    if (d_reader->GetImageIO()->GetNumberOfDimensions() == 3)
    {
      typedef Image<ScalarType, 3>               ImageType;
      typedef ImageFileReader <ImageType>        ImageReaderType;
      typedef ImageFileWriter <ImageType>        ImageWriterType;
      typedef RegionOfInterestImageFilter< ImageType, ImageType > ExtractorType;

      ImageReaderType::Pointer     reader          = ImageReaderType::New();
      ExtractorType::Pointer       extractor       = ExtractorType::New();
      ImageWriterType::Pointer     writer          = ImageWriterType::New();
      
      std::cout<<"reading input 3D file"<<std::endl;  
      reader->SetFileName (inputfile);
      reader->Update();
      std::cout<<"done."<<std::endl;
      ImageType::Pointer input = reader->GetOutput();
      
      std::cout<<"extracting region"<<std::endl;  
      ImageType::RegionType region = FindRegion<ImageType,TransformType>(input, transform, diameter);
      extractor->SetInput (input);
      extractor->SetRegionOfInterest (region);
      extractor->Update();
      std::cout<<"done."<<std::endl;
      ImageType::Pointer output = extractor->GetOutput();
      
      std::cout << "Writing image to : "<< outputfile << std::endl;
      writer->SetFileName (outputfile);
      writer->SetInput(output);
      writer->Update();
      std::cout<<"done."<<std::endl;
    }
    else if (d_reader->GetImageIO()->GetNumberOfDimensions() == 4)
    {
      typedef Image<ScalarType, 4>               ImageType;
      typedef ImageFileReader <ImageType>        ImageReaderType;
      typedef ImageFileWriter <ImageType>        ImageWriterType;
      typedef RegionOfInterestImageFilter< ImageType, ImageType > ExtractorType;

      ImageReaderType::Pointer     reader          = ImageReaderType::New();
      ExtractorType::Pointer       extractor       = ExtractorType::New();
      ImageWriterType::Pointer     writer          = ImageWriterType::New();
      
      std::cout<<"reading input 4D file"<<std::endl;  
      reader->SetFileName (inputfile);
      reader->Update();
      std::cout<<"done."<<std::endl;
      ImageType::Pointer input = reader->GetOutput();
      
      std::cout<<"extracting region"<<std::endl;  
      ImageType::RegionType region = FindRegion<ImageType,TransformType>(input, transform, diameter);
      extractor->SetInput (input);
      extractor->SetRegionOfInterest (region);
      extractor->Update();
      std::cout<<"done."<<std::endl;
      ImageType::Pointer output = extractor->GetOutput();
      
      std::cout << "Writing image to : "<< outputfile << std::endl;
      writer->SetFileName (outputfile);
      writer->SetInput(output);
      writer->Update();
      std::cout<<"done."<<std::endl;
    }
    
    return EXIT_SUCCESS;
  }

}
