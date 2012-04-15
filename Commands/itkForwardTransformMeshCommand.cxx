#include "itkForwardTransformMeshCommand.h"

#include "itkProlateSpheroidalTransformTensorMeshFilter.h"

#include "itkTensorImageIO.h"
#include "itkTensorMeshIO.h"
#include <itkImageFileReader.h>
#include <itkTensorImageToMeshFilter.h>
#include <itkLimitToAHAZoneImageFilter.h>

#include "itkWarpTensorMeshFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include <itksys/SystemTools.hxx>

#include <sstream>
#include <fstream>
#include <vector>
#include <fstream>
#include <iostream>

#include "GetPot.h"

namespace itk
{

  ForwardTransformMeshCommand::ForwardTransformMeshCommand()
  {
    m_ShortDescription = "Transform a tensor image or mesh in Prolate Spheroidal Coordinates";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input tensor image/mesh (default : input.vtk)]\n";    
    m_LongDescription +="-pr   [prolate transform used]\n";
    m_LongDescription +="-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription +="-o    [output vtk unstructured tensor mesh]\n";    
  }

  ForwardTransformMeshCommand::~ForwardTransformMeshCommand()
  {}

  int ForwardTransformMeshCommand::Execute (int narg, const char* arg[])
  {
    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout  << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char* inputfile                    = cl.follow("input.vtk",2,"-i","-I");
    const char* prolatefile                  = cl.follow("prolate.lms",2,"-pr","-PR");
    const char* displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char* inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char* outputfile                   = cl.follow("output.csv",2,"-o","-O");

    // typedefs
    typedef double                                                         ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>                           TensorImageIOType;
    typedef itk::TensorMeshIO <ScalarType, 3, 3>                           TensorMeshIOType;

    typedef TensorImageIOType::TensorImageType                             TensorImageType;
    typedef TensorImageType::PixelType                                     TensorType;

    typedef itk::Vector<ScalarType, 3>                                         DisplacementType;
    typedef itk::Image<DisplacementType, 3>                                DisplacementFieldType;
    typedef itk::ImageFileReader<DisplacementFieldType>                    DisplacementFileReaderType;

    typedef TensorMeshIOType::TensorMeshType                               MeshType;
    typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>     WarperType;
    typedef itk::TensorImageToMeshFilter<TensorType, 3>                    TensorImageToMeshFilterType;

    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>      TransformerType;
    typedef TransformerType::TransformType                                 TransformType;
    typedef TransformType::InputPointType                                  PointType;

    // instantiation
    DisplacementFileReaderType::Pointer    displacementreader1    = DisplacementFileReaderType::New();
    DisplacementFileReaderType::Pointer    displacementreader2    = DisplacementFileReaderType::New();
    TransformType::Pointer                 transform              = TransformType::New();
    TransformerType::Pointer               transformer            = TransformerType::New();
    WarperType::Pointer                    warper                 = WarperType::New();
    
    // read the input tensors and put tham into a vtkUnstructuredGrid
    // they come from a text file listing all files to read, either vtk or itk...
    MeshType::Pointer data = NULL;
    
    std::cout<<"reading input : "<<inputfile<<std::endl;
    std::string extension = itksys::SystemTools::GetFilenameLastExtension(inputfile).c_str();
    if (strcmp (extension.c_str(), ".vtk") == 0)
    {
      TensorMeshIOType::Pointer reader = TensorMeshIOType::New();
      reader->SetFileName (inputfile);
      reader->Read();
      data = reader->GetOutput();
    }
    else
    {
      TensorImageIOType::Pointer reader = TensorImageIOType::New();
      reader->SetFileName(inputfile);    
      reader->Read();
      TensorImageToMeshFilterType::Pointer imagetomesh = TensorImageToMeshFilterType::New();
      imagetomesh->SetInput (reader->GetOutput());
      imagetomesh->Update();
      data = imagetomesh->GetOutput();
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
  
    std::cout<<"reading transform "<<prolatefile<<"..."<<std::endl;
    itk::TransformFactory<TransformType>::RegisterTransform ();
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
  
    transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    TransformType::Pointer transform_inverse = TransformType::New();
    transform->GetInverse(transform_inverse);
    std::cout << " Done." << std::endl;
    
    std::cout<<"Warping..."<<std::endl;
    warper->SetInput (data);
    warper->SetDisplacementField (displacementfield);
    warper->SetInverseDisplacementField (inversedisplacementfield);
    warper->Update();
    
    std::cout<<"Transforming..."<<std::endl;
    transformer->SetInput (warper->GetOutput());
    transformer->SetTransform (transform);
    transformer->Update();
    
    MeshType::Pointer output = transformer->GetOutput();

    output->DisconnectPipeline();

    std::cout<<"Writing..."<<std::endl;
    TensorMeshIOType::Pointer writer = TensorMeshIOType::New();
    writer->SetFileName (outputfile);
    writer->SetInput (output);
    writer->Write();
    std::cout<<" Done."<<std::endl;

    return EXIT_SUCCESS;
  }
  
}
