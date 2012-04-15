/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkExtrapolateTensorField.cxx 1 2010-05-21 14:00:33Z nt08 $
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
#include "itkCreateProlateAtlasCommand.h"

#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkTensorImageIO.h"
#include "itkTensorMeshIO.h"
#include "itkProlateSpheroidalTransformTensorMeshFilter.h"
#include "itkWarpTensorMeshFilter.h"
#include <itkTensorImageToMeshFilter.h>
#include <itkTensorMeshToImageFilter.h>


#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include "GetPot.h"
#include <itkGaussianInterpolationTensorMeshFilter2.h>

template<typename MeshType, typename ImageType>
typename MeshType::Pointer DomainToMesh (typename ImageType::Pointer image)
{

  typename MeshType::Pointer mesh = MeshType::New();
  
  itk::ImageRegionIterator<ImageType>   itIn(image, image->GetLargestPossibleRegion());
  typename MeshType::PointsContainer::Pointer    points = MeshType::PointsContainer::New();
  typename MeshType::PointDataContainer::Pointer data   = MeshType::PointDataContainer::New();
  mesh->SetPoints (points);
  mesh->SetPointData (data);
    
  typename MeshType::PointType x;
  unsigned int counter = 0;
    
  while(!itIn.IsAtEnd())
  {
    image->TransformIndexToPhysicalPoint (itIn.GetIndex(), x);
    if (itIn.Get() > vcl_numeric_limits<typename ImageType::PixelType>::epsilon())
      mesh->SetPoint (counter++, x);
    ++itIn;
  }
    
  return mesh;
}

template<typename MeshType, typename ImageType>
void UnNormalizeFirstComponent(typename MeshType::Pointer in, double* range)
{

  typename MeshType::PointsContainer::Pointer points = in->GetPoints();
  typename MeshType::PointsContainer::Iterator it_points = points->Begin();

  double min = range[0];
  double max = range[1];
  
  std::cout<<"1st component unnormalization between : "<<min<<" and "<<max<<std::endl;
  
  it_points = points->Begin();
  unsigned long counter = 0;
  
  while(it_points != points->End())
  {
    typename MeshType::PointType p = it_points.Value();
    p[0] = min + (max - min) * (p[0]);
    in->SetPoint (counter, p);
    ++it_points;
    counter++;
  }
}

namespace itk
{

  CreateProlateAtlasCommand::CreateProlateAtlasCommand()
  {
    m_ShortDescription = "Create an tensor atlas from a set of tensors in prolate coordinates ";
    m_LongDescription += m_ShortDescription;
    m_LongDescription += "\nThe key parameter is the kernel sizes";
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-i    [input tensor unstructured grid (default : input.vtk)]\n";
    m_LongDescription += "-d    [input domaiin (default : domain.mha)]\n";
    m_LongDescription += "-pr   [prolate transform used (default : prolate.pr)]\n";
    m_LongDescription += "-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription += "-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription += "-u     [kernel size list file to use]\n";
    m_LongDescription += "-o    [output tensor atlas\n";
  }
  
  CreateProlateAtlasCommand::~CreateProlateAtlasCommand()
  {}

  int CreateProlateAtlasCommand::Execute (int narg, const char* arg[])
  {

    typedef double                                                         ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>                           TensorImageIOType;
    typedef itk::TensorMeshIO <ScalarType, 3, 3>                           TensorMeshIOType;
    
    typedef TensorImageIOType::TensorImageType                             TensorImageType;
    typedef itk::Image<ScalarType,3>                                       ImageType;
    typedef TensorImageType::PixelType                                     TensorType;
    typedef TensorMeshIOType::TensorMeshType                               MeshType;

    typedef itk::Vector<ScalarType, 3>                                     DisplacementType;
    typedef itk::Image<DisplacementType, 3>                                DisplacementFieldType;
    typedef itk::ImageFileReader<DisplacementFieldType>                    DisplacementFileReaderType;
    typedef itk::ImageFileWriter<DisplacementFieldType>                    DisplacementFileWriterType;
    typedef itk::ImageFileReader<ImageType>                                ImageFileReaderType;
    typedef itk::ImageFileWriter<ImageType>                                ImageFileWriterType;

    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>      TransformerType;
    typedef TransformerType::TransformType                                 TransformType;
    
    typedef MeshType::PointType                     PointType;
    typedef itk::ImageRegionIterator<ImageType>     ImageIteratorType;

    typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType> WarperType;
    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>  TransformerType;

    typedef itk::GaussianInterpolationTensorMeshFilter2<MeshType>  InterpolatorType;
    typedef itk::TensorMeshToImageFilter<TensorType, 3> TensorMeshToImageFilterType;

    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
  
    const char*  inputfile                    = cl.follow("NoFile",2,"-i","-I");
    const char*  domainfile                    = cl.follow("NoFile",2,"-d","-D");
    const char*  prolatefile                  = cl.follow("NoFile",2,"-pr","-PR");
    const char*  displacementfieldfile        = cl.follow("NoFile",2,"-f1","-F1");
    const char*  inversedisplacementfieldfile = cl.follow("NoFile",2,"-f2","-F2");
    const char*  outputfile                    = cl.follow("NoFile",2,"-o","-O");
    const char* kernelfile    = cl.follow("nofile", 2,"-u","-U");
  
    std::cout << "Processing atlas creation: " << std::endl;
    std::cout << std::flush;
  
    // read the input tensors and put tham into a vtkUnstructuredGrid
  
    std::cout << "Reading input tensors: " << inputfile << std::endl;
    TensorMeshIOType::Pointer reader = TensorMeshIOType::New();
    reader->SetFileName (inputfile);
    reader->Read();
    MeshType::Pointer data = reader->GetOutput();

    ImageFileReaderType::Pointer dreader = ImageFileReaderType::New();
    dreader->SetFileName (domainfile);
    dreader->Update();
    ImageType::Pointer domainimage = dreader->GetOutput();
  
  
    std::cout<<"reading transformation"<<std::endl;
    itk::TransformFactory<TransformType>::RegisterTransform ();
    TransformType::Pointer transform = TransformType::New();  
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
    transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
  
    DisplacementFileReaderType::Pointer    displacementreader1 = DisplacementFileReaderType::New();
    DisplacementFileReaderType::Pointer    displacementreader2 = DisplacementFileReaderType::New();

    // read the displacement field images
    DisplacementFieldType::Pointer displacementfield = NULL;
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


  
  
    std::cout<<"reading kernel list : "<<kernelfile<<std::endl;
    std::ifstream inputliststream (kernelfile);
    if(inputliststream.fail())
    {
      std::cerr << "Unable to open file: " << kernelfile << std::endl;
      std::exit (EXIT_FAILURE);
    }
    unsigned int NumberOfKernels = 0;
    inputliststream >> NumberOfKernels;
    std::cout<<"encountered : "<<NumberOfKernels<<std::endl;
  
    std::string sline = "";
    itksys::SystemTools::GetLineFromStream(inputliststream, sline);
    std::vector<double*> kernellist;
    for (unsigned int N=0; N<NumberOfKernels; N++)
    {
      std::string line = "";
      itksys::SystemTools::GetLineFromStream(inputliststream, line);
      itksys_ios::istringstream parse ( line );
      double* kernel = new double[3];
      for (unsigned int i=0; i<3; i++)
	parse >> kernel[i];
    
      kernellist.push_back (kernel);
    }

    std::cout<<"found "<<NumberOfKernels<<" kernels :"<<std::endl;
    double* alphas = new double[3*kernellist.size()];
    for (unsigned int N=0; N<NumberOfKernels; N++)
    {
      std::cout<<N<<": "
	       <<kernellist[N][0]<<" : "
	       <<kernellist[N][1]<<" : "
	       <<kernellist[N][2]<<std::endl;
      alphas[3*N+0] = 1.0 * kernellist[N][0];
      alphas[3*N+1] = 1.0 * kernellist[N][1];
      alphas[3*N+2] = 1.0 * kernellist[N][2];
    
    }
  

    WarperType::Pointer             WarperIn = WarperType::New();
    TransformerType::Pointer SwitcherIn = TransformerType::New();
    WarperIn->SetDisplacementField (displacementfield);
    WarperIn->SetInverseDisplacementField (inversedisplacementfield);
    SwitcherIn->SetTransform (transform);

    MeshType::Pointer domaincartesian = DomainToMesh<MeshType,ImageType> (domainimage);

    WarperIn->SetInput (domaincartesian);
    WarperIn->Update();
    SwitcherIn->SetInput (WarperIn->GetOutput());
    SwitcherIn->Update();
    MeshType::Pointer domain = SwitcherIn->GetOutput();

    MeshType::PointsContainer::Pointer points = domain->GetPoints();
    MeshType::PointsContainer::Iterator it_points = points->Begin();

    double max = vcl_numeric_limits<double>::min();
    double min = vcl_numeric_limits<double>::max();
    double d = 0.0;
  
    while(it_points != points->End())
    {
      d = it_points.Value()[0];    
      if (d < min) min = d;
      if (d > max) max = d;
      ++it_points;
    }
    double range[2] = {min, max};

    UnNormalizeFirstComponent<MeshType,ImageType> (data, range);

    TensorMeshIOType::Pointer writert = TensorMeshIOType::New();
    writert->SetInput (data);
    writert->SetFileName("data.vtk");
    writert->Write();
    writert->SetInput (domain);
    writert->SetFileName("domain.vtk");
    writert->Write();

    domain->DisconnectPipeline();
    data->DisconnectPipeline();
  
    InterpolatorType::Pointer m_Interpolator = InterpolatorType::New();

    m_Interpolator->SetInput (0, domain);
    m_Interpolator->SetInput (1, data);
    m_Interpolator->SetUsePiWorkAround(1);
    m_Interpolator->SetAlpha (alphas);
  
    InterpolatorType::LimiterType::Pointer limiter = m_Interpolator->GetLimiter();
    limiter->CanineDivisionsOff();
    limiter->SetTransform (transform);
    limiter->SetDisplacementField (displacementfield);
    limiter->SetInverseDisplacementField (inversedisplacementfield);
    limiter->SetAHASegmentationType (InterpolatorType::LimiterType::AHA_1_ZONE);
    limiter->CalculateZones();  

    m_Interpolator->Update();
  
    MeshType::Pointer output = m_Interpolator->GetOutput();
    output->DisconnectPipeline();

  
    writert->SetInput (output);
    writert->SetFileName("output.vtk");
    writert->Write();

    WarperType::Pointer             WarperOut = WarperType::New();
    TransformerType::Pointer SwitcherOut = TransformerType::New();
    WarperOut->SetDisplacementField (inversedisplacementfield);
    WarperOut->SetInverseDisplacementField (displacementfield);
    TransformType::Pointer backtransform = TransformType::New();
    transform->GetInverse (backtransform);
    SwitcherOut->SetTransform (backtransform);

    SwitcherOut->SetInput (data);
    SwitcherOut->Update();
  
    WarperOut->SetInput (SwitcherOut->GetOutput());
    WarperOut->Update();

    std::cout << "Writing the output... "<<outputfile<< " " << std::flush;
    std::string extension = itksys::SystemTools::GetFilenameLastExtension(outputfile).c_str();
    if (strcmp (extension.c_str(), ".vtk") == 0)
    {
      TensorMeshIOType::Pointer writer = TensorMeshIOType::New();
      writer->SetFileName (outputfile);
      writer->SetInput (WarperOut->GetOutput());
      writer->Write();
    }
    else
    {
      TensorMeshToImageFilterType::Pointer meshtoimage = TensorMeshToImageFilterType::New();
      meshtoimage->SetInput (WarperOut->GetOutput());
      meshtoimage->SetDomain (domainimage);
      meshtoimage->Update();
      TensorImageIOType::Pointer writer = TensorImageIOType::New();
      writer->SetFileName(outputfile);
      writer->SetInput (meshtoimage->GetOutput());
      writer->Write();
    }
    std::cout<<" Done."<<std::endl;

    SwitcherOut->SetInput (data);
    SwitcherOut->Update();
  
    WarperOut->SetInput (SwitcherOut->GetOutput());
    WarperOut->Update();
    TensorMeshIOType::Pointer writer = TensorMeshIOType::New();
    writer->SetFileName ("rawdata.vtk");
    writer->SetInput (WarperOut->GetOutput());
    writer->Write();

    return EXIT_SUCCESS;
  }

}
