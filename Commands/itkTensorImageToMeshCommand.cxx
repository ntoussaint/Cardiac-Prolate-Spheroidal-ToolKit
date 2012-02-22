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
#include "itkTensorImageToMeshCommand.h"

#include <itkImageFileReader.h>

#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>


#include <itkTensorImageToMeshFilter.h>
#include <itkTensor.h>
#include <itkTensorImageIO.h>
#include <itkTensorMeshIO.h>

#include "GetPot.h"


namespace itk
{

  TensorImageToMeshCommand::TensorImageToMeshCommand()
  {
    m_ShortDescription = "Converts a tensor image into a tensor mesh structure";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription +="-i    [input  image]\n";
    m_LongDescription +="-o    [output mesh]\n";
  }

  TensorImageToMeshCommand::~TensorImageToMeshCommand()
  {}

  int TensorImageToMeshCommand::Execute (int narg, const char* arg[])
  {

    typedef double                                      ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>        TensorImageIOType;
    typedef itk::TensorMeshIO<ScalarType, 3, 3>         TensorMeshIOType;
    typedef TensorImageIOType::TensorType               TensorType;
    typedef itk::TensorImageToMeshFilter<TensorType, 3> FilterType;
    typedef FilterType::MeshType                        MeshType;
    typedef FilterType::ImageType                       ImageType;
  
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char*  inputfile  = cl.follow("input.vtk",2,"-i","-I");
    const char*  outputfile = cl.follow("output.mha",2,"-o","-O");
  
    std::cout << "inputfile: \t\t\t" << inputfile << std::endl;
    std::cout << "outputfile: \t\t\t" <<outputfile << std::endl;

    std::cout << "Reading input tensors: " << inputfile <<"... "<< std::flush;  
    TensorImageIOType::Pointer reader = TensorImageIOType::New();
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
    ImageType::Pointer input = reader->GetOutput();
  
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput (reader->GetOutput());

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
    TensorMeshIOType::Pointer writer = TensorMeshIOType::New();
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
