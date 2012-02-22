/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkCreateSyntheticTensorMap.cxx 1 2010-05-21 14:00:33Z nt08 $
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
#include "itkCreateSyntheticCardiacTensorMapCommand.h"


/**
 * \file itkCreateSyntheticCardiacTensorMapCommand.cxx
 * 
 * \seealso ProlateSpheroidalTransformTensorMeshFilter WarpTensorMeshFilter ProlateSpheroidalTransform
 * \author Nicolas Toussaint
 */


#include "itkTensorImageIO.h"
#include "itkTensorMeshIO.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTensorImageToMeshFilter.h>
#include <itkTensorMeshToImageFilter.h>

#include "itkProlateSpheroidalTransformTensorMeshFilter.h"
#include "itkWarpTensorMeshFilter.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

#include <itksys/SystemTools.hxx>

#include <vtkUnstructuredGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellType.h>

#include <sstream>
#include <fstream>
#include <vector>
#include <fstream>
#include <iostream>

#include "GetPot.h"

// typedefs
typedef double                                                         ScalarType;
typedef itk::TensorImageIO<ScalarType, 3, 3>                           TensorImageIOType;
typedef itk::TensorMeshIO <ScalarType, 3, 3>                           TensorMeshIOType;

typedef TensorImageIOType::TensorImageType                             TensorImageType;
typedef itk::Image<ScalarType,3>                                       ImageType;
typedef TensorImageType::PixelType                                     TensorType;
typedef itk::ImageFileReader<ImageType>                                ImageFileReaderType;
typedef itk::ImageFileWriter<ImageType>                                ImageFileWriterType;
typedef TensorImageType::PixelType                                     TensorType;  
typedef itk::ImageRegionIterator<ImageType>                            ImageIteratorType;
typedef itk::ImageRegionIterator<TensorImageType>                      TensorIteratorType;
typedef itk::GradientImageFilter<ImageType>                            GradientImageFilterType;
typedef GradientImageFilterType::OutputPixelType                       CovariantVectorType;
typedef itk::Vector<double, 3>                                          VectorType;
typedef itk::Image<VectorType, 3>                                      VectorImageType;
typedef itk::Image<CovariantVectorType, 3>                             GradientImageType;
typedef itk::Matrix<ScalarType, 3, 3>                                  MatrixType;
typedef itk::LinearInterpolateImageFunction<ImageType, double>         InterpolatorType;
typedef itk::Vector<double, 3>                                         DisplacementType;
typedef itk::Image<DisplacementType, 3>                                DisplacementFieldType;
typedef itk::ImageFileReader<DisplacementFieldType>                    DisplacementFileReaderType;
typedef itk::ImageFileWriter<DisplacementFieldType>                    DisplacementFileWriterType;

typedef TensorMeshIOType::TensorMeshType                               MeshType;
typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>      TransformerType;
typedef itk::ProlateSpheroidalTransform<ScalarType>                    TransformType;
typedef TransformType::InputPointType                                  PointType;
typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>     WarperType;
typedef itk::TensorImageToMeshFilter<TensorType, 3>                    TensorImageToMeshFilterType;
typedef itk::TensorMeshToImageFilter<TensorType, 3>                    TensorMeshToImageFilterType;



MeshType::Pointer DomainToMesh2 (ImageType::Pointer image)
{

  MeshType::Pointer mesh = MeshType::New();
    
  itk::ImageRegionIterator<ImageType>   itIn(image, image->GetLargestPossibleRegion());
  MeshType::PointsContainer::Pointer    points = MeshType::PointsContainer::New();
  MeshType::PointDataContainer::Pointer data   = MeshType::PointDataContainer::New();
  mesh->SetPoints (points);
  mesh->SetPointData (data);
    
  MeshType::PointType x;
  unsigned int counter = 0;
    
  while(!itIn.IsAtEnd())
  {
    image->TransformIndexToPhysicalPoint (itIn.GetIndex(), x);
    if ( itIn.Get() > vcl_numeric_limits<float>::epsilon() )
    {
      TensorType T (0.0);
      mesh->SetPoint (counter, x);
      mesh->SetPointData (counter, T);
      counter++;
    }
      
    ++itIn;
  }
    
  return mesh;
}


void GetXi1Range(const MeshType::Pointer mesh, double* range)
{

  const MeshType::PointsContainer::Pointer points = mesh->GetPoints();
  MeshType::PointsContainer::ConstIterator it_points = points->Begin();

  double min = vcl_numeric_limits<double>::max();
  double max = vcl_numeric_limits<double>::min();
  double d = 0.0;
  
  while(it_points != points->End())
  {
    d = it_points.Value()[0]; 
    if (d < min) min = d;
    if (d > max) max = d;
    ++it_points;
  }

  range[0] = min;
  range[1] = max;
}

TensorType GetProlateTensor (PointType xi, double* range, double l_helix, double l_sheet, double e2)
{
  double xmin=range[0];
  double xmax=range[1];
  double amin=-l_helix * vnl_math::pi / 180;
  double amax= l_helix * vnl_math::pi / 180;
  
  double x = (xi[0] - xmin) / (xmax - xmin);
  double a = - (amax-amin) * x - amin;
  double b = std::sin ( 5.0 * xi[2] ) * l_sheet * ( 1.0 + std::cos (x * 2.0 * vnl_math::pi) ) * vnl_math::pi / 360;

  if (e2 < 0.21) e2 = 0.21;
  if (e2 > 0.99) e2 = 0.99;
  
  TensorType::MatrixType v, e, r, rota, rotb;
  v.Fill (0.0); rota.Fill (0.0); rotb.Fill (0.0);
  
  v[2][0] = 1; v[0][1] = 1; v[1][2] = 1;
  e[0][0] = 1.0; e[1][1] = e2; e[2][2] = 0.20;

  r = v * e * TensorType::MatrixType (v.GetTranspose());
  
  rota[0][0] =  1.0;
  rota[1][1] =   std::cos (-a);
  rota[2][1] = - std::sin (-a);
  rota[1][2] =   std::sin (-a);
  rota[2][2] =   std::cos (-a);  

  rotb[0][0] =   std::cos (b);
  rotb[1][0] = - std::sin (b);
  rotb[0][1] =   std::sin (b);
  rotb[1][1] =   std::cos (b);
  rotb[2][2] =   1.0;
  
  TensorType T;
  T.SetVnlMatrix ( r.GetVnlMatrix() );
  
  // std::cout<<"x = \n"<<x<<std::endl;
  // std::cout<<"a = \n"<<a<<std::endl;
  // std::cout<<"b = \n"<<b<<std::endl;
  // std::cout<<"v = \n"<<v<<std::endl;
  // std::cout<<"r = \n"<<r<<std::endl;
  // std::cout<<"rota = \n"<<rota<<std::endl;
  // std::cout<<"rotb = \n"<<rotb<<std::endl;
  // std::cout<<"cp : "<<T.GetCp()<<std::endl;
  // getchar();  
  
  
  T = T.ApplyMatrix (rotb);
  T = T.ApplyMatrix (rota);

  return T;
}


namespace itk
{

  CreateSyntheticCardiacTensorMapCommand::CreateSyntheticCardiacTensorMapCommand()
  {
    m_ShortDescription = "Create a synthetic tensor map out of some cardiac LV fibre structure information";
    m_LongDescription += m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-i    [input domain (default : domain.mha)]\n";
    m_LongDescription += "-pr   [prolate transform]\n";
    m_LongDescription += "-o    [output tensor map]\n";
    m_LongDescription += "-helix [helix angle max (default 60]\n";
    m_LongDescription += "-sheet [sheet angle max (default 60]\n";
    m_LongDescription += "-p [planar coefficient (0.2 < p < 0.9)]\n";;
  }

  CreateSyntheticCardiacTensorMapCommand::~CreateSyntheticCardiacTensorMapCommand()
  {}
  
  int CreateSyntheticCardiacTensorMapCommand::Execute (int narg, const char* arg[])
  {
    
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char* inputfile   = cl.follow("domain.mha",2,"-i","-I");
    const char* prolatefile = cl.follow("prolate.lms",2,"-pr","-PR");
    const char* outputfile  = cl.follow("output.mha",2,"-o","-O");
    const double l_helix    = cl.follow(60, 2,"-helix","-HELIX");
    const double l_sheet    = cl.follow(60, 2,"-sheet","-SHEET");
    const double p          = cl.follow(0.7, 2,"-p","-P");
  
    std::cout << "Processing bandwidth extraction with following arguments: " << std::endl;
    std::cout << "inputfile: " << inputfile << std::endl;
    std::cout << "prolatefile: " << prolatefile << std::endl;
    std::cout << "output: " << outputfile << std::endl;
    std::cout << "l_helix: " << l_helix<< std::endl;
    std::cout << "l_sheet: " << l_sheet<< std::endl;
    std::cout << "p: " << p<< std::endl;
  
    std::cout << std::flush;
  
    std::cout<<"reading input : "<<inputfile<<std::endl;
    ImageFileReaderType::Pointer domainreader = ImageFileReaderType::New();
    domainreader->SetFileName(inputfile);
    try
    {
      domainreader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << " Done." << std::endl;
    ImageType::Pointer domain = domainreader->GetOutput();  
    std::cout<<" Done."<<std::endl;

    MeshType::Pointer tensormesh = DomainToMesh2 (domain);
  
    itk::TransformFactory<TransformType>::RegisterTransform ();
  
    // instantiation
    TransformType::Pointer   transform   = TransformType::New();
    TransformerType::Pointer transformer = TransformerType::New();
  
    std::cout<<"reading"<<std::endl;
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
  
    transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    TransformType::Pointer transform_inverse = TransformType::New();
    transform->GetInverse(transform_inverse);
  
    transformer->SetInput (tensormesh);
    transformer->SetTransform (transform);
    transformer->Update();

    MeshType::Pointer prolatemesh = transformer->GetOutput();
    prolatemesh->DisconnectPipeline();
  
    double range[2];

    GetXi1Range (prolatemesh, range);
  
    MeshType::PointsContainer::Iterator it_points = prolatemesh->GetPoints()->Begin();
    MeshType::PointDataContainer::Iterator it_data   = prolatemesh->GetPointData()->Begin();
  
    while(it_points != prolatemesh->GetPoints()->End())
    {
      PointType xi = it_points.Value();

      TensorType T = GetProlateTensor (xi, range, l_helix, l_sheet, p);
      it_data.Value() = T;
    
      ++it_points;
      ++it_data;
    }

  
    TransformerType::Pointer backtransformer = TransformerType::New();
    TransformType::Pointer backtransform = TransformType::New();
    transform->GetInverse (backtransform);
  
    backtransformer->SetInput (prolatemesh);
    backtransformer->SetTransform (backtransform);
    backtransformer->Update();

    MeshType::Pointer output = backtransformer->GetOutput();
    output->DisconnectPipeline();
  
    std::cout << "Writing the output... "<<outputfile<< " " << std::flush;
    std::string extension = itksys::SystemTools::GetFilenameLastExtension(outputfile).c_str();
    if (strcmp (extension.c_str(), ".vtk") == 0)
    {
      TensorMeshIOType::Pointer writer = TensorMeshIOType::New();
      writer->SetFileName (outputfile);
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
      writer->SetFileName(outputfile);
      writer->SetInput (meshtoimage->GetOutput());
      writer->Write();
    }
    std::cout<<" Done."<<std::endl;

  
    return EXIT_SUCCESS;
  }

}
