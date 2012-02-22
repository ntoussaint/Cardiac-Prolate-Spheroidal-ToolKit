/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkExtractComponentDistribution.cxx 1 2010-05-21 14:00:33Z nt08 $
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
#include <itkWarpTensorMeshCommand.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

#include <itksys/SystemTools.hxx>

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>

#include "itkTensorImageIO.h"
#include "itkProlateSpheroidalTransformTensorMeshFilter.h"
#include "itkWarpTensorMeshFilter.h"

#include <sstream>
#include <fstream>
#include <vector>
#include <fstream>
#include <iostream>

#include "GetPot.h"

namespace itk
{

  WarpTensorMeshCommand::WarpTensorMeshCommand()
  {
    
    m_ShortDescription = "Warp";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription +="-i    [input tensor mesh (default : input.vtk)]\n";
    m_LongDescription +="-f1   [displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription +="-o    [output tensor mesh]\n";
  }


  WarpTensorMeshCommand::~WarpTensorMeshCommand()
  {}

  int WarpTensorMeshCommand::Execute (int narg, const char* arg[])
  {
  
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
  
    const char* inputfile                    = cl.follow("input.vtk",2,"-i","-I");
    const char* displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char* outputfile                   = cl.follow("output.vtk",2,"-o","-O");
  
    std::cout << "Processing bandwidth extraction with following arguments: " << std::endl;
    std::cout << "inputfile: " << inputfile << std::endl;
    std::cout << "displacementfieldfile: " << displacementfieldfile << std::endl;
    std::cout << "inversedisplacementfieldfile: " << inversedisplacementfieldfile << std::endl;
    std::cout << "output: " << outputfile << std::endl;
    std::cout << std::flush;

  
    // typedefs
    typedef double                                                         ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>                           TensorIOType;
    typedef TensorIOType::TensorImageType                                  TensorImageType;
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

    typedef itk::DefaultStaticMeshTraits<TensorType, 3, 3, double, double, TensorType> MeshTraits;
    typedef itk::Mesh<TensorType, 3, MeshTraits>                                       MeshType;
    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>                  CoordinateSwitcherType;
    typedef itk::ProlateSpheroidalTransform<ScalarType>                                TransformType;
    typedef TransformType::InputPointType                                              PointType;
    typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>                 WarperType;


  
    // instantiation
    DisplacementFileReaderType::Pointer    displacementreader1    = DisplacementFileReaderType::New();
    DisplacementFileReaderType::Pointer    displacementreader2    = DisplacementFileReaderType::New();
    WarperType::Pointer                    warperdata             = WarperType::New();
  
    // read the input tensors and put tham into a vtkUnstructuredGrid
    // they come from a text file listing all files to read, either vtk or itk...  
 
    std::cout<<"reading input : "<<inputfile<<std::endl;
    MeshType::Pointer Data = MeshType::New();
    MeshType::PointsContainer::Pointer DataPoints = MeshType::PointsContainer::New();
    MeshType::PointDataContainer::Pointer  DataPixels = MeshType::PointDataContainer::New();  

    vtkDataSetReader* reader = vtkDataSetReader::New();
    reader->SetFileName (inputfile);
    reader->Update();

    unsigned long NumberOfDataPoints = vtkPointSet::SafeDownCast (reader->GetOutput())->GetPoints()->GetNumberOfPoints();
    DataPoints->Reserve (NumberOfDataPoints);
    DataPixels->Reserve (NumberOfDataPoints);

    Data->SetPoints (DataPoints);
    Data->SetPointData (DataPixels);
    unsigned long counter = 0;
  
    for (unsigned long i=0; i<NumberOfDataPoints; i++)
    {
      double* tensor = NULL;
      double* point  = vtkPointSet::SafeDownCast (reader->GetOutput())->GetPoint (i);
    
      TensorType T;
      if (vtkPointSet::SafeDownCast (reader->GetOutput())->GetPointData()->GetTensors())
      {
	tensor = vtkPointSet::SafeDownCast (reader->GetOutput())->GetPointData()->GetTensors()->GetTuple (i);
      
	T[0] = tensor[0];
	T[1] = tensor[1];
	T[2] = tensor[4];
	T[3] = tensor[2];
	T[4] = tensor[5];
	T[5] = tensor[8];
      }
    
      PointType x;
      x[0] = point[0];
      x[1] = point[1];
      x[2] = point[2];

      Data->SetPoint (counter, x);
      Data->SetPointData (counter, T);
      counter++;
    }
    std::cout<<" Done."<<std::endl;

    // read the displacement field images
    DisplacementFieldType::Pointer displacementfield = NULL;
    DisplacementFieldType::Pointer inversedisplacementfield = NULL;

    std::cout << "Reading forward field: " << displacementfieldfile << std::flush;
    displacementreader1->SetFileName(displacementfieldfile);
    try
    {
      displacementreader1->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << " Done." << std::endl;
  
    displacementfield = displacementreader1->GetOutput();

    std::cout << "Reading backward field: " << inversedisplacementfieldfile << std::flush;
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
    std::cout << " Done." << std::endl;

    inversedisplacementfield = displacementreader2->GetOutput();

    warperdata->SetInput (Data);
    warperdata->SetDisplacementField (displacementfield);
    warperdata->SetInverseDisplacementField (inversedisplacementfield);
    warperdata->Update();

    vtkPoints* outputpoints = vtkPoints::New();
    unsigned long numberofpoints = warperdata->GetOutput()->GetNumberOfPoints();
    PointType point;
  
    for (unsigned long i=0; i<numberofpoints; i++)
    {
      warperdata->GetOutput()->GetPoint (i, &point);
      outputpoints->InsertNextPoint (point.GetDataPointer());
    }

    vtkDataSet* output = NULL;
    vtkDoubleArray* data = vtkDoubleArray::New();
    data->SetNumberOfComponents (9);
    data->SetNumberOfTuples (numberofpoints);

    TensorType t (0.0);
    if (vtkUnstructuredGrid::SafeDownCast (reader->GetOutput()))
    {
      vtkUnstructuredGrid* u_output  = vtkUnstructuredGrid::New();
      u_output->SetPoints (outputpoints);
      if (vtkUnstructuredGrid::SafeDownCast (reader->GetOutput())->GetCells())
	u_output->SetCells (VTK_LINE, vtkUnstructuredGrid::SafeDownCast (reader->GetOutput())->GetCells());

      for (unsigned int i=0; i<numberofpoints; i++)
      {
      
	warperdata->GetOutput()->GetPointData (i, &t);
	double vals[9];
	vals[0] = t[0];
	vals[1] = t[1];
	vals[2] = t[3];
	vals[3] = t[1];
	vals[4] = t[2];
	vals[5] = t[4];
	vals[6] = t[3];
	vals[7] = t[4];
	vals[8] = t[5];
      
	data->SetTuple (i, vals);
      }
      u_output->GetPointData()->SetTensors (data);
    
      output = u_output;
    }
    else if ((vtkPolyData::SafeDownCast (reader->GetOutput())))
    {
      vtkPolyData* p_output = vtkPolyData::New();
      p_output->SetPoints (outputpoints);
      if (vtkPolyData::SafeDownCast (reader->GetOutput())->GetPolys())
	p_output->SetPolys (vtkPolyData::SafeDownCast (reader->GetOutput())->GetPolys());
    
      for (unsigned int i=0; i<numberofpoints; i++)
      {
	warperdata->GetOutput()->GetPointData (i, &t);
	double vals[9];
	vals[0] = t[0];
	vals[1] = t[1];
	vals[2] = t[3];
	vals[3] = t[1];
	vals[4] = t[2];
	vals[5] = t[4];
	vals[6] = t[3];
	vals[7] = t[4];
	vals[8] = t[5];
      
	data->SetTuple (i, vals);
      }
      p_output->GetPointData()->SetTensors (data);
    
      output = p_output;
    }

  
    vtkDataSetWriter* writer = vtkDataSetWriter::New();
    writer->SetInput (output);
    writer->SetFileName (outputfile);
    writer->Update();

    std::cout<<" done."<<std::endl;

    reader->Delete();
    output->Delete();
    writer->Delete();
    outputpoints->Delete();
  
    return EXIT_SUCCESS;
  }
}
