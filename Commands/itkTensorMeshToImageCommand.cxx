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
#include "itkTensorMeshToImageCommand.h"

#include <itkImageFileReader.h>

#include <itkTensorMeshToImageFilter.h>
#include <itkTensor.h>
#include <itkTensorImageIO.h>
#include <itkTensorMeshIO.h>

#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include "GetPot.h"

namespace itk
{

  TensorMeshToImageCommand::TensorMeshToImageCommand()
  {
    m_ShortDescription = "\nConverts a tensor mesh structure into a tensor image\n\n";
    m_LongDescription += m_ShortDescription;
    m_LongDescription = "Usage:\n";
    
    m_LongDescription += "-i    [input  mesh]";
    m_LongDescription += "-d    [domain image]";
    m_LongDescription += "-o    [output image]";
  }

  TensorMeshToImageCommand::~TensorMeshToImageCommand()
  {}

  int TensorMeshToImageCommand::Execute (int narg, const char* arg[])
  {
    
    typedef double                                      ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>        TensorImageIOType;
    typedef itk::TensorMeshIO<ScalarType, 3, 3>         TensorMeshIOType;
    typedef TensorImageIOType::TensorType               TensorType;
    typedef itk::TensorMeshToImageFilter<TensorType, 3> FilterType;
    typedef FilterType::MeshType                        MeshType;
    typedef FilterType::DomainImageType                 DomainImageType;
    typedef itk::ImageFileReader<DomainImageType>       ImageReaderType;
    
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char*  inputfile  = cl.follow("input.vtk",2,"-i","-I");
    const char*  domainfile = cl.follow("domain.mha",2,"-d","-D");
    const char*  outputfile = cl.follow("output.mha",2,"-o","-O");
  
    std::cout << "inputfile: \t\t\t" << inputfile << std::endl;
    std::cout << "domainfile: \t\t\t" << domainfile << std::endl;
    std::cout << "outputfile: \t\t\t" <<outputfile << std::endl;

    std::cout << "Reading input tensors: " << inputfile <<"... "<< std::flush;  
    TensorMeshIOType::Pointer reader = TensorMeshIOType::New();
    reader->SetFileName(inputfile);
  
    ImageReaderType::Pointer dreader = ImageReaderType::New();
    dreader->SetFileName(domainfile);
  
    try
    {
      reader->Read();
      dreader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << " Done." << std::endl;
    MeshType::Pointer input = reader->GetOutput();
  
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput (reader->GetOutput());
    filter->SetDomain (dreader->GetOutput());

    try
    {
      filter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
  
    std::cout << "writing output to "<<outputfile<<" ... "<< std::flush;
    TensorImageIOType::Pointer writer = TensorImageIOType::New();
    writer->SetFileName(outputfile);
    writer->SetInput (filter->GetOutput());
  
    try
    {
      writer->Write();
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
