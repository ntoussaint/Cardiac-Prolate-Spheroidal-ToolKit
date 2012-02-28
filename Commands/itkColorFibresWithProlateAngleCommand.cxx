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
#include "itkColorFibresWithProlateAngleCommand.h"

#include "itkTensorImageIO.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "itkProlateSpheroidalTransformTensorMeshFilter.h"
#include "itkWarpTensorMeshFilter.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

#include <itksys/SystemTools.hxx>

#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellType.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>

#include <sstream>
#include <fstream>
#include <vector>
#include "GetPot.h"
#include <fstream>
#include <iostream>

void CalculateVectorFromAngleHSV (double angle, double* range, double* anglevector)
{
  angle = (angle - range[0])/(range[1] - range[0]);

  angle = std::max (0.0, angle);
  angle = std::min (1.0, angle);
  //fix me !
  // we wanna use the inverse of the jet colorbar
  angle = 1.0 - angle;
  
  if (angle >= 0 && angle < 0.2)
  {
    anglevector[0] = 1 - (angle - 0.0) / 0.2;
    anglevector[1] = 0;
    anglevector[2] = 1;
  }
  else if (angle >= 0.2 && angle < 0.4)
  {
    anglevector[0] = 0;
    anglevector[1] = (angle - 0.2) / 0.2;
    anglevector[2] = 1;
  }
  else if (angle >= 0.4 && angle < 0.6)
  {
    anglevector[0] = 0;
    anglevector[1] = 1;
    anglevector[2] = 1 - (angle - 0.4) / 0.2;
  }
  else if (angle >= 0.6 && angle < 0.8)
  {
    anglevector[0] = (angle - 0.6) / 0.2;
    anglevector[1] = 1;
    anglevector[2] = 0;
  }
  else if (angle >= 0.8 && angle <= 1.0)
  {
    anglevector[0] = 1;
    anglevector[1] = 1 - (angle - 0.8) / 0.2;
    anglevector[2] = 0;
  }
}


void CalculateVectorFromAngleJet (double angle, double* range, double* anglevector)
{
  angle = (angle - range[0])/(range[1] - range[0]);

  angle = std::max (0.0, angle);
  angle = std::min (1.0, angle);
  //fix me !
  // we wanna use the inverse of the jet colorbar
  //angle = 1.0 - angle;
  
  if (angle >= 0 && angle < 0.125)
  {
    anglevector[0] = 0.5 + (angle - 0.0) / 0.25;
    anglevector[1] = 0;
    anglevector[2] = 0;
  }
  else if (angle >= 0.125 && angle < 0.375)
  {
    anglevector[0] = 1;
    anglevector[1] = 0.0 + (angle - 0.125) / 0.25;
    anglevector[2] = 0;
  }
  else if (angle >= 0.375 && angle < 0.625)
  {
    anglevector[0] = 1.0 - (angle - 0.375) / 0.25;
    anglevector[1] = 1;
    anglevector[2] = 0.0 + (angle - 0.375) / 0.25;
  }
  else if (angle >= 0.625 && angle < 0.875)
  {
    anglevector[0] = 0;
    anglevector[1] = 1.0 - (angle - 0.625) / 0.25;
    anglevector[2] = 1;
  }
  else if (angle >= 0.875 && angle <= 1.0)
  {
    anglevector[0] = 0;
    anglevector[1] = 0;
    anglevector[2] = 1.0 - (angle - 0.875) / 0.25;
  }
}

namespace itk
{

  ColorFibresWithProlateAngleCommand::ColorFibresWithProlateAngleCommand()
  {
    m_ShortDescription = "Colorify a fibre field with Prolate angular information";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i  [input fibre field]\n";
    m_LongDescription +="-t  [prolate transform]\n";
    m_LongDescription +="-f1 [forward displacement field]\n";
    m_LongDescription +="-f2 [forward displacement field]\n";    
    m_LongDescription +="-g  [global color (default: 0)]\n";    
    m_LongDescription +="-rmin  [desired HSV range min (default: -60)]\n";
    m_LongDescription +="-rmax  [desired HSV range max (default: +60)]\n";
    m_LongDescription +="-o  [output colored fibre field]\n";
    m_LongDescription +="-t  [type of the color-coding (default: helix)]\n";
    m_LongDescription +="available types:\n";
    m_LongDescription +="\t helix [helix angle (default)] \n";
    m_LongDescription +="\t transverse [transverse angle] \n";
    m_LongDescription +="\t fa [Fractional Anisotropy]\n";
    m_LongDescription +="\t cl-cp-cs [Linear / Planar and Spherical coefficients]\n";
    m_LongDescription +="\t sheet [sheet angle]\n";
    
  }

  ColorFibresWithProlateAngleCommand::~ColorFibresWithProlateAngleCommand()
  {}

  int ColorFibresWithProlateAngleCommand::Execute (int narg, const char* arg[])
  {

    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char* inputfile                    = cl.follow("input.vtk",2,"-i","-I");
    const char* prolatefile                  = cl.follow("prolate.lms",2,"-pr","-PR");
    const char* displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char* outputfile                   = cl.follow("output.fib",2,"-o","-O");
    const double rangemin                    = cl.follow(-60.0,2,"-rmin","-RMIN");
    const double rangemax                    = cl.follow(60.0,2,"-rmax","-RMAX");
    const bool globalcolor                   = cl.follow(0,2,"-g","-G");
    const char* typestring                   = cl.follow("helix",2,"-t","-T");

    unsigned int type = 0;

    std::cout << "Processing colorification with: " << std::endl;
    std::cout << "inputfile: " << inputfile << std::endl;
    std::cout << "prolatefile: " << prolatefile << std::endl;
    std::cout << "displacementfieldfile: " << displacementfieldfile << std::endl;
    std::cout << "inversedisplacementfieldfile: " << inversedisplacementfieldfile << std::endl;
    std::cout << "output: " << outputfile << std::endl;
    std::cout << "range: [" << rangemin <<":"<<rangemax<< std::endl;
  
    std::cout << std::flush;  

    
    if      (std::strcmp (typestring,"helix") == 0 )      type = 0;
    else if (std::strcmp (typestring,"transverse") == 0 ) type = 1;
    else if (std::strcmp (typestring,"fa") == 0 )         type = 2;
    else if (std::strcmp (typestring,"cl-cp-cs") == 0 )   type = 3;
    else if (std::strcmp (typestring,"sheet") == 0 )      type = 4;
    
    std::cout<<"computing the color-coding of ";
    switch(type)
    {
	case 0:
	default:
	  std::cout<<"helix angle";
	  break;
	case 1:
	  std::cout<<"transverse angle";
	  break;
	case 2:
	  std::cout<<"fractional anisotropy";
	  break;
	case 3:
	  std::cout<<"linear/planar/spherical coefficients";
	  break;
	case 4:
	  std::cout<<"sheet angle";
	  break;
    }

    std::cout<<"..."<<std::endl;
    
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
    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>                        CoordinateSwitcherType;
    typedef itk::ProlateSpheroidalTransform<ScalarType>                                            TransformType;
    typedef TransformType::InputPointType                                              PointType;
    typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>     WarperType;

  
    // instantiation
    DisplacementFileReaderType::Pointer    displacementreader1    = DisplacementFileReaderType::New();
    DisplacementFileReaderType::Pointer    displacementreader2    = DisplacementFileReaderType::New();
    DisplacementFileWriterType::Pointer    displacementwriter     = DisplacementFileWriterType::New();
    TransformType::Pointer                 transform              = TransformType::New();
    CoordinateSwitcherType::Pointer        coordinateswitcherdata = CoordinateSwitcherType::New();
    CoordinateSwitcherType::Pointer        coordinateswitcherref  = CoordinateSwitcherType::New();
    WarperType::Pointer                    warperdata             = WarperType::New();
    WarperType::Pointer                    warperref              = WarperType::New();
  
    // read the input tensors and put tham into a vtkUnstructuredGrid
    // they come from a text file listing all files to read, either vtk or itk...  
 
    std::cout<<"reading input : "<<inputfile<<std::endl;
    MeshType::Pointer Data = MeshType::New();
    MeshType::PointsContainer::Pointer DataPoints = MeshType::PointsContainer::New();
    MeshType::PointDataContainer::Pointer  DataPixels = MeshType::PointDataContainer::New();  

    vtkDataSetReader* reader = vtkDataSetReader::New();
    reader->SetFileName (inputfile);
    reader->Update();
    vtkPolyData* input = vtkPolyData::SafeDownCast (reader->GetOutput());
  
    unsigned long NumberOfDataPoints = input->GetPoints()->GetNumberOfPoints();
    DataPoints->Reserve (NumberOfDataPoints);
    DataPixels->Reserve (NumberOfDataPoints);
    
    Data->SetPoints (DataPoints);
    Data->SetPointData (DataPixels);
    unsigned long counter = 0;
    vtkFloatArray* tensors = vtkFloatArray::SafeDownCast (input->GetPointData()->GetArray ("Tensors"));
  
    for (unsigned long i=0; i<NumberOfDataPoints; i++)
    {
      double* tensor = tensors->GetTuple (i);
      double* point  = input->GetPoint (i);
    
      TensorType T;
      T[0] = tensor[0];
      T[1] = tensor[1];
      T[2] = tensor[2];
      T[3] = tensor[3];
      T[4] = tensor[4];
      T[5] = tensor[5];

      PointType x;
      x[0] = point[0];
      x[1] = point[1];
      x[2] = point[2];
      if (T.GetTrace() > 0.1)
      {
	Data->SetPoint (counter, x);
	Data->SetPointData (counter, T);
	counter++;
      }
    
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
  
    std::cout<<"registering"<<std::endl;
    itk::TransformFactory<TransformType>::RegisterTransform ();
  
    std::cout<<"reading"<<std::endl;
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
  
    transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    TransformType::Pointer transform_inverse = TransformType::New();
    transform->GetInverse(transform_inverse);
    
    coordinateswitcherdata->SetInput (warperdata->GetOutput());
    coordinateswitcherdata->SetTransform (transform);
    coordinateswitcherdata->Update();
  
    coordinateswitcherref->SetInput (coordinateswitcherdata->GetOutput());
    coordinateswitcherref->SetTransform (transform_inverse);
    coordinateswitcherref->Update();
  
    warperref->SetInput (coordinateswitcherref->GetOutput());
    warperref->SetDisplacementField (inversedisplacementfield);
    warperref->SetInverseDisplacementField (displacementfield);
    warperref->Update();

    vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
    colors->SetNumberOfComponents (3);
    colors->SetNumberOfTuples (input->GetNumberOfPoints());
    colors->FillComponent (0, (unsigned char)(255.0));
    colors->FillComponent (1, (unsigned char)(0.0));
    colors->FillComponent (2, (unsigned char)(0.0));
  
    double range[2]={rangemin, rangemax};
  
    unsigned long numberofpoints = Data->GetNumberOfPoints();
    PointType mainpointprolate;
    VectorType mainvectorprolate, lastvectorprolate;
    
    double helix, sheet, transverse, fa;
    MeshType::PixelType prolatetensor (0.0);

    for (unsigned int i=0; i<3; i++)
    {
      mainpointprolate[i] = 0.0;
    }
  
    std::ostringstream os;
    std::cout<<"how many points are we dealing with ?? "<<numberofpoints<<std::endl;
    MeshType::Pointer ProlateData   = coordinateswitcherdata->GetOutput();
  
    vtkCellArray* lines = input->GetLines();
    lines->InitTraversal();
    vtkIdType npt  = 0;
    vtkIdType *pto = 0;

    
    vtkIdType test = lines->GetNextCell (npt, pto);
    while( test!=0 )
    {

      double meandata = 0.0;
      double data = 0.0;
      
      for( int fi=0; fi<npt; fi++)
      {
	int i = pto[fi];
      
	ProlateData->GetPoint (i, &mainpointprolate);
	ProlateData->GetPointData (i, &prolatetensor);
      
	mainvectorprolate = prolatetensor.GetEigenvector (2);
	lastvectorprolate = prolatetensor.GetEigenvector (0);
	mainvectorprolate.Normalize();
	lastvectorprolate.Normalize();
	
	switch (type)
	{
	    case 0:
	    default:
	      helix = mainvectorprolate[1];
	      if (mainvectorprolate[2] < 0) helix = -helix;
	      helix = std::asin(helix) * 180.0 / vnl_math::pi;
	      data = helix;
	      break;
	    case 1:
	      transverse = mainvectorprolate[0];
	      if (mainvectorprolate[2] < 0) transverse = -transverse;
	      transverse = std::asin(transverse) * 180.0 / vnl_math::pi;
	      data = transverse;
	      break;
	    case 2:
	      fa = prolatetensor.GetFA();
	      data = fa;
	      break;
	    case 3:
	      data = prolatetensor.GetCl();
	      break;
	    case 4:
	      sheet = lastvectorprolate[0];
	      sheet = std::asin(sheet) * 180.0 / vnl_math::pi;
	      data = sheet;
	      break;
	}
	
	meandata = meandata + data;
	
	// for color
	if (!globalcolor)
	{
	  double anglevector[3];
	  CalculateVectorFromAngleJet (data, range, anglevector);
	  
	  for (unsigned int k=0; k<3; k++)
	  {
	    double c = fabs (anglevector[k])*255.0;
	    colors->SetComponent( i, k, (unsigned char)( c>255.0?255.0:c ) );
	  }
	}
      
      }

      meandata = meandata / (double)(npt);
    
      // for color
      if (globalcolor)
      {
	double anglevector[3];
	CalculateVectorFromAngleJet (meandata, range, anglevector);
      
	for( int fi=0; fi<npt; fi++)
	{
	  int i = pto[fi];
	
	  for (unsigned int k=0; k<3; k++)
	  {
	    double c = fabs (anglevector[k])*255.0;
	    colors->SetComponent( i, k, (unsigned char)( c>255.0?255.0:c ) );
	  }
	}
      }
     
      test = lines->GetNextCell (npt, pto);
    }

    std::cout << "Writing the fibres... "<<outputfile<< " " << std::endl;
    vtkPolyData* helixanglefibres = vtkPolyData::New();
    helixanglefibres->DeepCopy (input);
    helixanglefibres->GetPointData()->SetScalars (colors);
    
    vtkDataSetWriter* writer = vtkDataSetWriter::New();
    writer->SetInput (helixanglefibres);
    writer->SetFileName (outputfile);
    //writer->SetFileTypeToBinary();
    writer->Update();
  
    reader->Delete();
    writer->Delete();
    helixanglefibres->Delete();
    colors->Delete();
    
    return EXIT_SUCCESS;
  }

}
