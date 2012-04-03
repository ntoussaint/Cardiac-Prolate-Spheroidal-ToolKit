#include "itkExtractProlateInformationCommand.h"

#include "itkProlateSpheroidalTransformTensorMeshFilter.h"

#include "itkTensorImageIO.h"
#include "itkTensorMeshIO.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTensorImageToMeshFilter.h>
#include <itkLimitToAHAZoneImageFilter.h>

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

namespace itk
{

  ExtractProlateInformationCommand::ExtractProlateInformationCommand()
  {
    m_ShortDescription = "Extract meaningful information in prolate spheroidal coordinates";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input tensor image (default : input.vtk)]\n";    
    m_LongDescription +="-pr   [prolate transform used]\n";
    m_LongDescription +="-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription +="-o    [output csv file where cost function values are stored (default: output.csv)]\n";
    
    m_LongDescription +="-t    [type of the output (default: helix)]\n";
    m_LongDescription +="available types:\n";
    m_LongDescription +="\t helix [helix angle] \n";
    m_LongDescription +="\t transverse [transverse angle] \n";
    m_LongDescription +="\t fa [Fractional Anisotropy]\n";
    m_LongDescription +="\t cl-cp-cs [Linear / Planar and Spherical coefficients]\n";
    m_LongDescription +="\t sheet [sheet angle]\n";
    m_LongDescription +="\t error [position error accumulation]\n";
    m_LongDescription +="\t zone [AHA zone of point]\n";
    m_LongDescription +="\t deltaangle [absolute angular error accumulation] \n";
    m_LongDescription +="\t h1 [scale factor in xi1 direction] \n";
    m_LongDescription +="\t h2 [scale factor in xi2 direction] \n";
    m_LongDescription +="\t h3 [scale factor in xi3 direction] \n";
  }

  ExtractProlateInformationCommand::~ExtractProlateInformationCommand()
  {}

  int ExtractProlateInformationCommand::Execute (int narg, const char* arg[])
  {
    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout  << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }

    
    const char* inputfile                    = cl.follow("input.vtk",2,"-i","-I");
    const char* prolatefile                  = cl.follow("prolate.lms",2,"-pr","-PR");
    // const char* domainfile                   = cl.follow("nofile",2,"-d","-D");
    const char* displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char* outputfile                   = cl.follow("output.csv",2,"-o","-O");
    const char* typestring                   = cl.follow("helix",2,"-t","-T");

    unsigned int type = 0;

    if      (std::strcmp (typestring,"helix") == 0 )      type = 0;
    else if (std::strcmp (typestring,"transverse") == 0 ) type = 1;
    else if (std::strcmp (typestring,"fa") == 0 )         type = 2;
    else if (std::strcmp (typestring,"cl-cp-cs") == 0 )   type = 3;
    else if (std::strcmp (typestring,"sheet") == 0 )      type = 4;
    else if (std::strcmp (typestring,"error") == 0 )      type = 5;
    else if (std::strcmp (typestring,"zone") == 0 )       type = 6;
    else if (std::strcmp (typestring,"deltaangle") == 0 ) type = 7;
    else if (std::strcmp (typestring,"h1") == 0 )         type = 8;
    else if (std::strcmp (typestring,"h2") == 0 )         type = 9;
    else if (std::strcmp (typestring,"h3") == 0 )         type = 10;

    
    std::cout<<"computing the extraction of ";
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
	case 5:
	  std::cout<<"position error accumulation";
	  break;
	case 6:
	  std::cout<<"AHA zone";
	  break;
	case 7:
	  std::cout<<"angular error";
	  break;
	case 8:
	  std::cout<<"scale factor h1";
	  break;
	case 9:
	  std::cout<<"scale factor h2";
	  break;
	case 10:
	  std::cout<<"scale factor h3";
	  break;
    }

    std::cout<<"..."<<std::endl;
    
    
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
    typedef itk::Vector<ScalarType, 3>                                          VectorType;
    typedef itk::Image<VectorType, 3>                                      VectorImageType;
    typedef itk::Image<CovariantVectorType, 3>                             GradientImageType;
    typedef itk::Matrix<ScalarType, 3, 3>                                  MatrixType;
    typedef itk::LinearInterpolateImageFunction<ImageType, ScalarType>         InterpolatorType;
    typedef itk::Vector<ScalarType, 3>                                         DisplacementType;
    typedef itk::Image<DisplacementType, 3>                                DisplacementFieldType;
    typedef itk::ImageFileReader<DisplacementFieldType>                    DisplacementFileReaderType;
    typedef itk::ImageFileWriter<DisplacementFieldType>                    DisplacementFileWriterType;

    typedef TensorMeshIOType::TensorMeshType                               MeshType;
    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>      CoordinateSwitcherType;
    typedef CoordinateSwitcherType::TransformType                         TransformType;
    typedef TransformType::InputPointType                                  PointType;
    typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>     WarperType;
    typedef itk::TensorImageToMeshFilter<TensorType, 3>                    TensorImageToMeshFilterType;
    typedef itk::LimitToAHAZoneImageFilter<TensorImageType> AHALimiterType;

    // instantiation
    DisplacementFileReaderType::Pointer    displacementreader1    = DisplacementFileReaderType::New();
    DisplacementFileReaderType::Pointer    displacementreader2    = DisplacementFileReaderType::New();
    DisplacementFileWriterType::Pointer    displacementwriter     = DisplacementFileWriterType::New();
    TransformType::Pointer                 transform              = TransformType::New();
    CoordinateSwitcherType::Pointer        coordinateswitcherdata = CoordinateSwitcherType::New();
    CoordinateSwitcherType::Pointer        coordinateswitcherref  = CoordinateSwitcherType::New();
    WarperType::Pointer                    warperdata             = WarperType::New();
    WarperType::Pointer                    warperref              = WarperType::New();

    TensorImageType::Pointer tensorimage = 0;
    
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

      tensorimage = reader->GetOutput();
      
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

    AHALimiterType::Pointer zonelimiter = AHALimiterType::New();
    zonelimiter->SetInput (tensorimage);
    zonelimiter->SetAHASegmentationType (AHALimiterType::AHA_17_ZONES);
    zonelimiter->SetTransform (transform);
    zonelimiter->SetInverseDisplacementField (inversedisplacementfield);
    zonelimiter->CanineDivisionsOff();
    
    try
    {
      zonelimiter->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }

    
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
  
    std::cout << "Writing csv file in : " << outputfile << std::endl;
    unsigned long numberofpoints = Data->GetNumberOfPoints();
    VectorType mainvector, lastvector, error, initialvector;
    VectorType mainvectorprolate, lastvectorprolate;
    PointType mainpoint; mainpoint[0] = mainpoint[1] = mainpoint[2] = 0.0;
    PointType mainpointprolate; mainpointprolate[0] = mainpointprolate[1] = mainpointprolate[2] = 0.0;
    PointType initialpoint; initialpoint[0] = initialpoint[1] = initialpoint[2] = 0.0;
    std::ostringstream os;
    std::cout<<"how many points are we dealing with ?? "<<numberofpoints<<std::endl;
    MeshType::Pointer ProlateData   = coordinateswitcherdata->GetOutput();
    MeshType::Pointer CartesianData = warperref->GetOutput();
    
    TensorType prolatetensor (0.0);
    TensorType cartesiantensor (0.0);
    TensorType initialtensor (0.0);

    double helix=0, sheet=0, transverse=0, fa=0, deltav = 0, scalefactor[3]={0,0,0};
    unsigned int zone=0;
    
    for (unsigned long i=0; i<numberofpoints; i++)
    {
      CartesianData->GetPointData (i, &cartesiantensor);
      CartesianData->GetPoint (i, &mainpoint);
      ProlateData->GetPointData (i, &prolatetensor);
      ProlateData->GetPoint (i, &mainpointprolate);
      Data->GetPoint (i, &initialpoint);
      Data->GetPointData (i, &initialtensor);
      
      mainvector = cartesiantensor.GetEigenvector (2);
      mainvectorprolate = prolatetensor.GetEigenvector (2);
      lastvector = cartesiantensor.GetEigenvector (0);
      lastvectorprolate = prolatetensor.GetEigenvector (0);
      initialvector = initialtensor.GetEigenvector (2);
      
      mainvector.Normalize();
      mainvectorprolate.Normalize();
      lastvector.Normalize();
      lastvectorprolate.Normalize();
      initialvector.Normalize();

      transform->EvaluateScaleFactors (mainpointprolate.GetDataPointer(), scalefactor);
      
      os << mainpointprolate[0] << " "
	 << mainpointprolate[1] << " "
	 << mainpointprolate[2] << " ";
      
      switch (type)
      {
	  case 0:
	  default:
	    helix = mainvectorprolate[1];
	    if (mainvectorprolate[2] < 0) helix = -helix;
	    helix = std::asin(helix) * 180.0 / vnl_math::pi;
	    os << helix;
	    break;
	  case 1:
	    transverse = mainvectorprolate[0];
	    if (mainvectorprolate[2] < 0) transverse = -transverse;
	    transverse = std::asin(transverse) * 180.0 / vnl_math::pi;
	    os << transverse;
	    break;
	  case 2:
	    fa = prolatetensor.GetFA();
	    os << fa;
	    break;
	  case 3:
	    os << prolatetensor.GetCl() << " ";
	    os << prolatetensor.GetCp() << " ";
	    os << prolatetensor.GetCs();
	    break;
	  case 4:
	    sheet = lastvectorprolate[0];
	    sheet = std::asin(sheet) * 180.0 / vnl_math::pi;
	    os << sheet;
	    break;
	  case 5:
	    error = mainpoint - initialpoint;
	    os << error.GetNorm();
	    break;
	  case 6:
	    zone = zonelimiter->InWhichZoneIsCartesianPoint (mainpoint);
	    os << zone;
	    break;
	  case 7:
	    deltav = std::abs (mainvector * initialvector);
	    os << std::acos (deltav) * 180 / vnl_math::pi;
	    break;
	  case 8:
	    os << scalefactor[0];
	    break;
	  case 9:
	    os << scalefactor[1];
	    break;
	  case 10:
	    os << scalefactor[2];
	    break;
      }
      
      os << std::endl;
      
      positions->SetTuple (i, mainpointprolate.GetDataPointer());
      projections->SetTuple (i, mainvectorprolate.GetDataPointer());
      mainvectors->SetTuple (i, mainvector.GetDataPointer());
      positionerrors->SetTuple (i,error.GetDataPointer());
      coordinateswitcherref->GetOutput()->GetPoint (i, &pt);
      ellipsoidpoints->SetPoint (i, pt[0], pt[1], pt[2]);
      pt = mainpoint;
      ventriclepoints->SetPoint (i, pt[0], pt[1], pt[2]);
      helixarray->SetTuple1 (i, helix);
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
  
    // vtkDataSetWriter* writer = vtkDataSetWriter::New();
    // writer->SetInput (outputgrid);
    // writer->SetFileName ("ellipsoiddada.vtk");
    // writer->SetFileTypeToBinary();
    // writer->Update();

    // outputgrid->SetPoints (ventriclepoints);
    // writer->SetInput (outputgrid);
    // writer->SetFileName ("ventricledata.vtk");
    // writer->SetFileTypeToBinary();
    // writer->Update();
    // writer->Delete();

    outputgrid->Delete();
    ellipsoidpoints->Delete();
    ventriclepoints->Delete();
  
    std::cout<<" done."<<std::endl;

    return EXIT_SUCCESS;
  }
  
}
