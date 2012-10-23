#include "itkForwardTransformImageCommand.h"

#include "itkProlateSpheroidalTransformTensorMeshFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

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

  ForwardTransformImageCommand::ForwardTransformImageCommand()
  {
    m_ShortDescription = "Transform a scalar image in a Prolate Spheroidal Coordinates box";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input scalar image (default : input.mha)]\n";    
    m_LongDescription +="-pr   [prolate transform used]\n";
    m_LongDescription +="-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription +="-o    [output image in prolate coordinates]\n";    
  }

  ForwardTransformImageCommand::~ForwardTransformImageCommand()
  {}

  int ForwardTransformImageCommand::Execute (int narg, const char* arg[])
  {
    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout  << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char* inputfile                    = cl.follow("input.mha",2,"-i","-I");
    const char* prolatefile                  = cl.follow("prolate.lms",2,"-pr","-PR");
    const char* displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char* outputfile                   = cl.follow("output.csv",2,"-o","-O");
    
    // define the future box boundaries :
    
    double PSS_box_bounds[3][2];
    
    PSS_box_bounds[0][0] = 0.3;
    PSS_box_bounds[0][1] = 0.6;
    
    PSS_box_bounds[1][0] = 0.0;
    PSS_box_bounds[1][1] = vnl_math::pi / 2.0;
    
    PSS_box_bounds[2][0] = 0.0;
    PSS_box_bounds[2][1] = 2.0 * vnl_math::pi;

    // define the future box size :

    unsigned int PSS_box_size[3];

    PSS_box_size[0] = 50;
    PSS_box_size[1] = 100;
    PSS_box_size[2] = 100;
    
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

    
    std::cout<<"reading input : "<<inputfile<<std::endl;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inputfile);    
    reader->Update();  
    ImageType::Pointer inputimage = reader->GetOutput();
    std::cout<<" Done."<<std::endl;

    // instantiation
    DisplacementFileReaderType::Pointer    displacementreader1    = DisplacementFileReaderType::New();
    DisplacementFileReaderType::Pointer    displacementreader2    = DisplacementFileReaderType::New();
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
    
    ImageType::DirectionType direction;
    direction.SetIdentity();
    ImageType::PointType origin;
    ImageType::SpacingType spacing;
    ImageType::RegionType region;
    ImageType::SizeType size;
    origin[0] = PSS_box_bounds[0][0];
    origin[1] = PSS_box_bounds[1][0];
    origin[2] = PSS_box_bounds[2][0];
    size[0] = PSS_box_size[0];
    size[1] = PSS_box_size[1];
    size[2] = PSS_box_size[2];
    for (unsigned int i=0; i<3; i++)
      spacing[i] = (double)(PSS_box_bounds[i][1] - PSS_box_bounds[i][0])/(double)(PSS_box_size[i] + 1);
    region.SetSize (size);
    ImageType::Pointer outputimage = ImageType::New();
    outputimage->SetRegions (region);
    outputimage->SetOrigin(origin);
    outputimage->SetSpacing(spacing);  
    outputimage->SetDirection(direction);
    outputimage->Allocate();
    outputimage->FillBuffer (0.0);

    // iterate and fill the image
    
    // for each point in the output image xi, do
    // p = phi * psi (xi)
    // ask input image the scalar value at p
    // put the value in output image
    
    itk::ImageRegionIterator<ImageType> itOut (outputimage, outputimage->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType> itIn  (inputimage,  inputimage->GetLargestPossibleRegion());
    ImageType::PointType xi, p1, p2;
    itk::ContinuousIndex<ScalarType, 3> index;

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetInputImage (inputimage);
    DisplacementInterpolatorType::Pointer displacementinterpolator = DisplacementInterpolatorType::New();
    displacementinterpolator->SetInputImage (inversedisplacementfield);
    
    while (!itOut.IsAtEnd())
    {
      ImageType::PixelType value = static_cast<ImageType::PixelType>(0.0);
      
      outputimage->TransformIndexToPhysicalPoint (itOut.GetIndex(), xi);
      p1 = transform_inverse->TransformPoint (xi); // psi-1 operator
      if (displacementinterpolator->IsInsideBuffer (p1))
      {
	p2 = p1 + displacementinterpolator->Evaluate (p1); // phi-1 operator

	if (interpolator->IsInsideBuffer (p2))
	  value = interpolator->Evaluate (p2); // ask the input image its value at p2 and put it in the box
      }
      itOut.Set (value);
      ++itOut;
    }

    // to trick the spacing so that the image is more viewable...
    spacing[0] *= 200;
    spacing[1] *= 90;
    spacing[2] *= 30;

    outputimage->SetSpacing (spacing);
    
	
    std::cout<<"Writing..."<<std::endl;
    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName (outputfile);
    writer->SetInput (outputimage);
    writer->Update();
    std::cout<<" Done."<<std::endl;

    return EXIT_SUCCESS;
  }
  
}
