#include "itkBackwardTransformImageCommand.h"

#include "itkProlateSpheroidalTransformTensorMeshFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include <itksys/SystemTools.hxx>
#include "itkVectorLinearInterpolateImageFunction.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <fstream>
#include <iostream>

#include "GetPot.h"

namespace itk
{

  BackwardTransformImageCommand::BackwardTransformImageCommand()
  {
    m_ShortDescription = "Transform Prolate Spheroidal Coordinates box into a scalar image";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input scalar image (default : input.mha)]\n";    
    m_LongDescription +="-pr   [prolate transform used]\n";
    m_LongDescription +="-f    [BACKWARD displacement field (default : forward.mha)]\n";
    m_LongDescription +="-o    [output image in prolate coordinates]\n";    
  }

  BackwardTransformImageCommand::~BackwardTransformImageCommand()
  {}

  int BackwardTransformImageCommand::Execute (int narg, const char* arg[])
  {
    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout  << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char* inputfile                    = cl.follow("input.mha",2,"-i","-I");
    const char* prolatefile                  = cl.follow("prolate.lms",2,"-pr","-PR");
    const char* displacementfieldfile        = cl.follow("backward.mha",2,"-f","-F");
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
    typedef itk::LinearInterpolateImageFunction<ImageType, ScalarType>  InterpolatorType;
    typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
    
    std::cout<<"reading input : "<<inputfile<<std::endl;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inputfile);    
    reader->Update();  
    ImageType::Pointer inputimage = reader->GetOutput();
    std::cout<<" Done."<<std::endl;
    
    // instantiation
    DisplacementFileReaderType::Pointer    displacementreader1    = DisplacementFileReaderType::New();
    // read the displacement field images
    DisplacementFieldType::Pointer displacementfield = NULL;

    std::cout << "Reading backward field: " << displacementfieldfile << std::flush;
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
    // we take the bounds of the displacement field,
    // they should englobe the LV.
    ImageType::DirectionType direction;
    ImageType::PointType origin;
    ImageType::SpacingType spacing;
    ImageType::RegionType region;
    ImageType::SizeType size;
    direction = displacementfield->GetDirection();
    origin = displacementfield->GetOrigin();
    spacing = displacementfield->GetSpacing();
    size = displacementfield->GetLargestPossibleRegion().GetSize();
    region.SetSize (size);
    ImageType::Pointer outputimage = ImageType::New();
    outputimage->SetRegions (region);
    outputimage->SetOrigin(origin);
    outputimage->SetSpacing(spacing);  
    outputimage->SetDirection(direction);
    outputimage->Allocate();
    outputimage->FillBuffer (0.0);

    // iterate and fill the image
    
    // for each point in the cartesian output image p, do
    // xi = psi * phi (p)
    // ask input image the scalar value at xi
    // put the value in output cartesian image at p
    
    itk::ImageRegionIterator<ImageType> itOut (outputimage, outputimage->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType> itIn  (inputimage,  inputimage->GetLargestPossibleRegion());
    ImageType::PointType xi, p1, p2;
    itk::ContinuousIndex<ScalarType, 3> index;

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetInputImage (inputimage);
    DisplacementInterpolatorType::Pointer displacementinterpolator = DisplacementInterpolatorType::New();
    displacementinterpolator->SetInputImage (displacementfield);
    
    while (!itOut.IsAtEnd())
    {
      ImageType::PixelType value = static_cast<ImageType::PixelType>(0.0);
      
      outputimage->TransformIndexToPhysicalPoint (itOut.GetIndex(), p1);
      if (displacementinterpolator->IsInsideBuffer (p1))
      {
	p2 = p1 + displacementinterpolator->Evaluate (p1); // phi operator	
      }
      xi = transform->TransformPoint (p2); // psi operator
      if (interpolator->IsInsideBuffer (xi))
	value = interpolator->Evaluate (xi); // ask the input image its value at p2 and put it in the box
      
      itOut.Set (value);
      ++itOut;
    }
    
    std::cout<<"Writing..."<<std::endl;
    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName (outputfile);
    writer->SetInput (outputimage);
    writer->Update();
    std::cout<<" Done."<<std::endl;

    return EXIT_SUCCESS;
  }
  
}
