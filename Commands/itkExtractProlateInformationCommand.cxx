#include "itkExtractProlateInformationCommand.h"
#include <iostream>
#include "GetPot.h"

#include "itkProlateSpheroidalTransformTensorMeshFilter.h"

#include "itkTensorImageIO.h"
#include "itkTensorMeshIO.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTensorImageToMeshFilter.h>

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
#include "GetPot.h"
#include <fstream>
#include <iostream>

namespace itk
{

  ExtractProlateInformationCommand::ExtractProlateInformationCommand()
  {
    m_ShortDescription = "Extract meaningful information in prolate spheroidal coordinates\n\n";
    m_LongDescription = "Usage:\n";
    m_LongDescription +="-i    [input tensor image (default : input.vtk)]\n";    
    m_LongDescription +="-pr   [prolate transform used]\n";
    m_LongDescription +="-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription +="-o    [output csv file where cost function values are stored (default: output.csv)]\n";
    m_LongDescription +="-s    [sheet-angle distribution instead of elevation (default 0)]\n";

    m_LongDescription += m_ShortDescription;
  }

  ExtractProlateInformationCommand::~ExtractProlateInformationCommand()
  {}

  int ExtractProlateInformationCommand::Execute (int narg, const char* arg[])
  {

    itk::Object::GlobalWarningDisplayOff();
    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << this->GetLongDescription() << std::endl;
      return -1;
    }

    
  const char* inputfile                    = cl.follow("input.vtk",2,"-i","-I");
  const char* prolatefile                  = cl.follow("prolate.lms",2,"-pr","-PR");
  // const char* domainfile                   = cl.follow("nofile",2,"-d","-D");
  const char* displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
  const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
  const char* outputfile                   = cl.follow("output.csv",2,"-o","-O");
  const bool  sheetangle                   = cl.follow(false,2,"-s","-S");
  
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
typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>      CoordinateSwitcherType;
 typedef CoordinateSwitcherType::TransformType                         TransformType;
typedef TransformType::InputPointType                                  PointType;
typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>     WarperType;
typedef itk::TensorImageToMeshFilter<TensorType, 3>                    TensorImageToMeshFilterType;



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
  std::string extension = itksys::SystemTools::GetFilenameLastExtension(inputfile).c_str();
  MeshType::Pointer Data = NULL;
  if (strcmp (extension.c_str(), ".vtk") == 0)
  {
    TensorMeshIOType::Pointer reader = TensorMeshIOType::New();
    reader->SetFileName (inputfile);
    reader->Read();
    Data = reader->GetOutput();
  }
  else
  {
    TensorImageIOType::Pointer reader = TensorImageIOType::New();
    reader->SetFileName(inputfile);    
    reader->Read();
    TensorImageToMeshFilterType::Pointer imagetomesh = TensorImageToMeshFilterType::New();
    imagetomesh->SetInput (reader->GetOutput());
    imagetomesh->Update();
    Data = imagetomesh->GetOutput();
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

  vtkDoubleArray* positions = vtkDoubleArray::New();
  positions->SetNumberOfComponents (3);
  positions->SetNumberOfTuples (Data->GetNumberOfPoints());
  positions->SetName ("PPS coord");

  vtkDoubleArray* positionerrors = vtkDoubleArray::New();
  positionerrors->SetNumberOfComponents (3);
  positionerrors->SetNumberOfTuples (Data->GetNumberOfPoints());
  positionerrors->SetName ("position errors");
  
  vtkDoubleArray* mainvectors = vtkDoubleArray::New();
  mainvectors->SetNumberOfComponents (3);
  mainvectors->SetNumberOfTuples (Data->GetNumberOfPoints());
  mainvectors->SetName ("mainvectors");

  vtkDoubleArray* projections = vtkDoubleArray::New();
  projections->SetNumberOfComponents (3);
  projections->SetNumberOfTuples (Data->GetNumberOfPoints());
  projections->SetName ("PPS proj");

  vtkDoubleArray* helixarray = vtkDoubleArray::New();
  helixarray->SetNumberOfComponents (1);
  helixarray->SetNumberOfTuples (Data->GetNumberOfPoints());
  helixarray->SetName ("helix angle");
  
  vtkPoints* ellipsoidpoints = vtkPoints::New();
  ellipsoidpoints->SetNumberOfPoints (Data->GetNumberOfPoints());
  vtkPoints* ventriclepoints = vtkPoints::New();
  ventriclepoints->SetNumberOfPoints (Data->GetNumberOfPoints());

  PointType pt; pt[0] = pt[1] = pt[2] = 0.0;
  double pos[3] = {0,0,0};
  double err[3] = {0,0,0};  
  double proj[3] = {0,0,0};
  double vec[3] = {0,0,0};
  
  std::cout << "Writing csv file in : " << outputfile << std::endl;
  unsigned long numberofpoints = Data->GetNumberOfPoints();
  VectorType mainvector;
  VectorType mainvectorprolate;
  PointType mainpoint; mainpoint[0] = mainpoint[1] = mainpoint[2] = 0.0;
  PointType mainpointprolate; mainpointprolate[0] = mainpointprolate[1] = mainpointprolate[2] = 0.0;
  PointType initialpoint; initialpoint[0] = initialpoint[1] = initialpoint[2] = 0.0;
  std::ostringstream os;
  std::cout<<"how many points are we dealing with ?? "<<numberofpoints<<std::endl;
  MeshType::Pointer ProlateData   = coordinateswitcherdata->GetOutput();
  MeshType::Pointer CartesianData = warperref->GetOutput();
  
  TensorType prolatetensor (0.0);
  TensorType cartesiantensor (0.0);
  
  for (unsigned long i=0; i<numberofpoints; i++)
  {
    ProlateData->GetPointData (i, &prolatetensor);
    CartesianData->GetPointData (i, &cartesiantensor);

    if (sheetangle)
    {
      mainvector = cartesiantensor.GetEigenvector (0);
      mainvectorprolate = prolatetensor.GetEigenvector (0);
    }
    else
    {  
      mainvector = cartesiantensor.GetEigenvector (2);
      mainvectorprolate = prolatetensor.GetEigenvector (2);
    }
    
    CartesianData->GetPoint (i, &mainpoint);
    ProlateData->GetPoint (i, &mainpointprolate);
    Data->GetPoint (i, &initialpoint);

    mainvector.Normalize();
    mainvectorprolate.Normalize();

    
    for (unsigned int j=0; j<3; j++)
    {
      pos[j] = mainpointprolate[j];
      vec[j] = mainvector[j];
      err[j] = mainpoint[j] - initialpoint[j];      
      // @warning The eigenvectors are assumed to be normalized.
      proj[j] = mainvectorprolate[j];
    }
    
#if 0
    os << pos[0] << " "
       << pos[1] << " "
       << pos[2] << " "
       << proj[0] << " "
       << proj[1] << " "
       << proj[2] << " "
       << std::endl;
#endif
    
#if 1
    // VectorType secondvector = prolatetensor.GetEigenvector (1);
    // secondvector.Normalize();
    // double anisotropy  = std::acos (secondvector[0]);    
    // os << anisotropy << " " << std::endl;
    double c1 = prolatetensor.GetCl();
    double cp = prolatetensor.GetCp();
    double cs = 1 - (c1 + cp);
    os << pos[0] << " "
       << pos[1] << " "
       << pos[2] << " ";
    os << c1 << " "
       << cp << " "
       << cs << std::endl;
#endif
    
    positions->SetTuple (i, pos);
    projections->SetTuple (i, proj);
    mainvectors->SetTuple (i, vec);
    positionerrors->SetTuple (i,err);
    coordinateswitcherref->GetOutput()->GetPoint (i, &pt);
    ellipsoidpoints->SetPoint (i, pt[0], pt[1], pt[2]);
    pt = mainpoint;
    ventriclepoints->SetPoint (i, pt[0], pt[1], pt[2]);

    double angle = proj[0];
    if (sheetangle)
    {
      if (proj[1])
	angle = -angle;
    }
    else
    {  
      if (proj[2])
	angle = -angle;
    }
    angle = std::asin(angle) * 180.0 / vnl_math::pi;
    helixarray->SetTuple1 (i, angle);
    
  }


  std::ofstream buffer (outputfile, ios::out);
  buffer << os.str().c_str();
  buffer.close();  

  std::cout<<" done."<<std::endl;

  std::cout<<" writing the output plain vtk grid."<<std::endl;

  vtkUnstructuredGrid* outputgrid = vtkUnstructuredGrid::New();
  outputgrid->SetPoints (ellipsoidpoints);
  outputgrid->GetPointData()->AddArray (positions);
  outputgrid->GetPointData()->AddArray (positionerrors);
  outputgrid->GetPointData()->AddArray (projections);
  outputgrid->GetPointData()->AddArray (helixarray);
  outputgrid->GetPointData()->SetVectors (mainvectors);
  
  positions->Delete();
  positionerrors->Delete();
  projections->Delete();
  mainvectors->Delete();
  
  vtkDataSetWriter* writer = vtkDataSetWriter::New();
  writer->SetInput (outputgrid);
  writer->SetFileName ("ellipsoiddada.vtk");
  writer->SetFileTypeToBinary();
  writer->Update();

  outputgrid->SetPoints (ventriclepoints);
  writer->SetInput (outputgrid);
  writer->SetFileName ("ventricledata.vtk");
  writer->SetFileTypeToBinary();
  writer->Update();
  writer->Delete();

  outputgrid->Delete();
  ellipsoidpoints->Delete();
  ventriclepoints->Delete();
  
  std::cout<<" done."<<std::endl;

  return EXIT_SUCCESS;
  }
  
}
