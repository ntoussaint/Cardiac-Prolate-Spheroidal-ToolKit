#include "itkCreateProlateCoordinateImageCommand.h"

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
    
  CreateProlateCoordinateImageCommand::CreateProlateCoordinateImageCommand()
  {
    m_ShortDescription = "Simply Create an image with 3-Component Prolate Coordinates from input domain";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input (tensor) image]\n";    
    m_LongDescription +="-pr   [prolate transform used]\n";
    m_LongDescription +="-f   [BACKWARD displacement field (default : backward.mha)]\n";
    m_LongDescription +="-o    [output 3-component image]\n";        
  }
    
  CreateProlateCoordinateImageCommand::~CreateProlateCoordinateImageCommand()
  {}
    
  int CreateProlateCoordinateImageCommand::Execute (int narg, const char* arg[])
  {   
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout  << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
      
    const char* inputfile                    = cl.follow("input.mha",2,"-i","-I");
    const char* prolatefile                  = cl.follow("prolate.tr",2,"-pr","-PR");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f","-F");
    const char* outputfile                   = cl.follow("output.csv",2,"-o","-O");
    
    // typedefs
    typedef double                                                         ScalarType;
    typedef itk::Image<ScalarType,3>                                       ImageType;
    typedef itk::ImageFileReader<ImageType>                                ImageReaderType;    
    typedef itk::ImageFileWriter<ImageType>                                ImageWriterType;    
    typedef itk::Vector<ScalarType, 3>                                     DisplacementType;
    typedef itk::Image<DisplacementType, 3>                                DisplacementFieldType;
    typedef itk::ImageFileReader<DisplacementFieldType>                    DisplacementFileReaderType;
    typedef itk::ProlateSpheroidalTransform<ScalarType>                    TransformType;
    typedef TransformType::InputPointType                                  PointType;
    typedef VectorLinearInterpolateImageFunction<DisplacementFieldType,ScalarType> DisplacementInterpolatorType;
    typedef itk::Image<ScalarType, 4>                                      Image4DType;
    typedef itk::ImageFileWriter<Image4DType>                              Image4DFileWriterType;
    typedef itk::ImageRegionConstIterator<ImageType>                       ImageIteratorType;
    typedef itk::ImageRegionIterator<Image4DType>                          Image4DIteratorType;
    
    std::cout<<"reading input : "<<inputfile<<std::endl;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inputfile);    
    reader->Update();  
    ImageType::Pointer inputimage = reader->GetOutput();
    std::cout<<" Done."<<std::endl;
      
    // instantiation
    DisplacementFileReaderType::Pointer    displacementreader2    = DisplacementFileReaderType::New();
    // read the displacement field images
    DisplacementFieldType::Pointer inversedisplacementfield = NULL;
      
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
    
    std::cout<<"reading transform "<<prolatefile<<"..."<<std::endl;
    itk::TransformFactory<TransformType>::RegisterTransform ();
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
      
    TransformType::Pointer transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    TransformType::Pointer transform_inverse = TransformType::New();
    transform->GetInverse(transform_inverse);
    std::cout << " Done." << std::endl;
    
    // create the output image geometry -> outputimage
    Image4DType::Pointer outputimage = Image4DType::New();
    Image4DType::RegionType region;
    Image4DType::SizeType size;
    Image4DType::SpacingType spacing;
    Image4DType::PointType origin;
    Image4DType::DirectionType direction;
    direction.Fill (0.0);
    
    for (unsigned int i=0; i<3; i++)
    {
      size[i] = inputimage->GetLargestPossibleRegion().GetSize()[i];
      spacing[i] = inputimage->GetSpacing()[i];
      origin[i] = inputimage->GetOrigin()[i];
      for (unsigned int j=0; j<3; j++)
	direction[i][j] = inputimage->GetDirection()[i][j];
    }
    size[3] = 3;
    spacing[3] = 1.0;
    origin[3] = 0.0;
    direction[3][3] = 1.0;
    region.SetSize (size);
    
    outputimage->SetRegions (region);
    outputimage->SetOrigin(origin);
    outputimage->SetSpacing(spacing);  
    outputimage->SetDirection(direction);
    outputimage->Allocate();
    outputimage->FillBuffer (0.0);
    
    itk::ImageRegionIterator<Image4DType> itOut (outputimage, outputimage->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType> itIn  (inputimage,  inputimage->GetLargestPossibleRegion());
    
    itk::ContinuousIndex<ScalarType, 3> continuousindex;
    Image4DType::IndexType outputindex;
    ImageType::IndexType index;
    
    PointType x;
    
    DisplacementInterpolatorType::Pointer displacementinterpolator = DisplacementInterpolatorType::New();
    displacementinterpolator->SetInputImage (inversedisplacementfield);
      
    while(!itIn.IsAtEnd())
    {
      if (!itIn.Get())
      {
	if(!itIn.IsAtEnd())
	  ++itIn;
	continue;
      }
      
      inputimage->TransformIndexToPhysicalPoint (itIn.GetIndex(), x);
      
      bool isinside = displacementinterpolator->IsInsideBuffer (x);
      if (!isinside)
      {
	if(!itIn.IsAtEnd())
	  ++itIn;
	continue;
      }
      inversedisplacementfield->TransformPhysicalPointToIndex (x, index);
      inversedisplacementfield->TransformPhysicalPointToContinuousIndex (x, continuousindex);
      DisplacementType d = displacementinterpolator->EvaluateAtContinuousIndex (continuousindex);
      
      x += d;
      
      PointType coordinates = transform->TransformPoint(x);
      
      for (unsigned int i=0; i<3; i++)
	outputindex[i] = itIn.GetIndex()[i];
      
      for (unsigned int i=0; i<3; i++)
      {
	outputindex[3] = i;
	itOut.SetIndex (outputindex);
	itOut.Set (coordinates[i]);
      }
      ++itIn;
    }
    Image4DFileWriterType::Pointer writer = Image4DFileWriterType::New();
    writer->SetInput (outputimage);
    writer->SetFileName (outputfile);
    writer->Update();
    
    
    return EXIT_SUCCESS;
  }
  
}
