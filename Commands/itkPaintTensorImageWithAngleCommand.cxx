#include "itkPaintTensorImageWithAngleCommand.h"

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

  PaintTensorImageWithAngleCommand::PaintTensorImageWithAngleCommand()
  {
    m_ShortDescription = "Paint meaningful information in prolate spheroidal coordinates";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input tensor image (default : input.vtk)]\n";    
    m_LongDescription +="-pr   [prolate transform used]\n";
    m_LongDescription +="-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription +="-o    [output mha file where function values are stored (default: output.mha)]\n";
    
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
    m_LongDescription +="\t l2 [L2 norm of the tensors] \n";
    m_LongDescription +="\t ta [transverse anisotropy] \n";
  }

  PaintTensorImageWithAngleCommand::~PaintTensorImageWithAngleCommand()
  {}

  int PaintTensorImageWithAngleCommand::Execute (int narg, const char* arg[])
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
    const char* outputfile                   = cl.follow("output.mha",2,"-o","-O");
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
    else if (std::strcmp (typestring,"l2") == 0 )         type = 11;
    else if (std::strcmp (typestring,"ta") == 0 )         type = 12;

    
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
	case 11:
	  std::cout<<"L2 norm";
	  break;
	case 12:
	  std::cout<<"transverse anisotropy";
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

    // read the input tensors and put tham into a vtkUnstructuredGrid
    // they come from a text file listing all files to read, either vtk or itk...  
    
    
    std::cout<<"reading input : "<<inputfile<<std::endl;
    std::string extension = itksys::SystemTools::GetFilenameLastExtension(inputfile).c_str();
    TensorImageType::Pointer tensorimage = 0;
    MeshType::Pointer Data = NULL;
    if (strcmp (extension.c_str(), ".mha") == 0)
    {
      TensorImageIOType::Pointer reader = TensorImageIOType::New();
      reader->SetFileName(inputfile);    
      reader->Read();

      tensorimage = reader->GetOutput();
      
      TensorImageToMeshFilterType::Pointer imagetomesh = TensorImageToMeshFilterType::New();
      imagetomesh->SetInput (reader->GetOutput());
      imagetomesh->Update();
      Data = imagetomesh->GetOutput();
    } else {
      std::cout<<"problem..."<<std::endl;
      std::exit (EXIT_FAILURE);
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
    MeshType::Pointer ProlateData   = coordinateswitcherdata->GetOutput();

    std::cout<<" evaluating the paint values."<<std::endl;
    // create the output image geometry -> outputimage
    // we take the bounds of the tensor field,
    ImageType::DirectionType direction;
    ImageType::PointType     origin;
    ImageType::SpacingType   spacing;
    ImageType::RegionType    region;
    ImageType::SizeType      size;
    direction = tensorimage->GetDirection();
    origin    = tensorimage->GetOrigin();
    spacing   = tensorimage->GetSpacing();
    size      = tensorimage->GetLargestPossibleRegion().GetSize();
    region.SetSize (size);
    ImageType::Pointer outputimage = ImageType::New();
    outputimage->SetRegions (region);
    outputimage->SetOrigin(origin);
    outputimage->SetSpacing(spacing);  
    outputimage->SetDirection(direction);
    outputimage->Allocate();
    outputimage->FillBuffer (0.0);
    
    VectorType mainvectorprolate;
    TensorType prolatetensor (0.0);
    double helix=0;
    unsigned int point_id = 0;
    itk::ImageRegionIterator<ImageType> itOut (outputimage, outputimage->GetLargestPossibleRegion());
    itk::ImageRegionIterator<TensorImageType> itIn  (tensorimage,  tensorimage->GetLargestPossibleRegion());
    
    while (!itOut.IsAtEnd())
    {
      if (itIn.Get().GetTrace() > vcl_numeric_limits<ScalarType>::epsilon())
      {
	ProlateData->GetPointData (point_id, &prolatetensor);
	mainvectorprolate = prolatetensor.GetEigenvector (2);
	mainvectorprolate.Normalize();
	switch (type)
	{
	    case 0:
	    default:
	      helix = mainvectorprolate[1];
	      if (mainvectorprolate[2] > 0) helix = -helix;
	      helix = std::asin(helix) * 180.0 / vnl_math::pi;
	      break;
	}
	
	itOut.Set (helix);
	point_id++;
      }
      
      ++itOut; ++itIn;
    }
    
    std::cout<<"Writing..."<<std::endl;
    ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
    writer->SetFileName (outputfile);
    writer->SetInput (outputimage);
    writer->Update();
    std::cout<<" Done."<<std::endl;

    return EXIT_SUCCESS;
  }
  
}
