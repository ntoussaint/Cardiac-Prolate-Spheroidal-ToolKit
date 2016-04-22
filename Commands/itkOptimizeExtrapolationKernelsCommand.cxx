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
#include "itkOptimizeExtrapolationKernelsCommand.h"

#include "itkTensorMeshImageHybridCostFunction3.h"
#include "itkNewUoaOptimizer.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkTensorImageIO.h"
#include "itkTensorMeshIO.h"
#include "itkTensorImageToMeshFilter.h"

#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include "GetPot.h"

namespace itk
{

  OptimizeExtrapolationKernelsCommand::OptimizeExtrapolationKernelsCommand()
  {
    
    m_ShortDescription = "Find optimal extrapolation Prol. Sph. kernel sizes to minimize discrepancies in a LS sense";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-p    [use prolate spheroid (default : 0.0)]\n";
    m_LongDescription += "-i    [input tensor unstructured grid (default : input.vtk)]\n";
    m_LongDescription += "-r    [reference tensors (default : reference.vtk)]\n";
    m_LongDescription += "-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription += "-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription += "-pr   [prolate transform used (default : prolate.pr)]\n";
    m_LongDescription += "-l    [lambda to use for smoothing ratio (default should be 0.01 but here set to 1.0)]\n";
    m_LongDescription += "-rhobegin [rho begin for the optimization (default is 0.03)]\n";
    m_LongDescription += "-rhoend   [rho end for the optimization (default is 0.005)]\n";
    m_LongDescription += "-u    [initial position for the optimizer\n";
    m_LongDescription += "-k    [interpolation kernel to use [0: Gaussian][1: B-Spline][2: Kaiser-Bessel](default: 0)]\n";
    m_LongDescription += "-o    [output csv file]\n\n";
    m_LongDescription += "CAUTION : it is important to create the outfile prior to execution\n";
    m_LongDescription += "          As the file is not created but the results are appended to it.\n";
  }


  OptimizeExtrapolationKernelsCommand::~OptimizeExtrapolationKernelsCommand()
  {}

  int OptimizeExtrapolationKernelsCommand::Execute (int narg, const char* arg[])
  {
  
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }

    typedef double                                  ScalarType;
    typedef itk::NewUoaOptimizer                    OptimizerType;
    typedef OptimizerType::ParameterVectorType      ParameterVectorType;
    typedef itk::TensorMeshImageHybridCostFunction3 CostFunctionType;
    typedef itk::Image<ScalarType,3>                ImageType;
    typedef CostFunctionType::TensorType            TensorType;
    typedef CostFunctionType::MeshType              MeshType;
    typedef CostFunctionType::DisplacementFieldType     DisplacementFieldType;
    typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFileReaderType;
  
    typedef CostFunctionType::TransformType         TransformType;
    typedef CostFunctionType::ParametersType        ParametersType;
    typedef CostFunctionType::ParametersValueType   ParametersValueType;
    typedef CostFunctionType::MeasureType           MeasureType;
    typedef MeshType::PointType                     PointType;
    typedef itk::ImageFileReader<ImageType>         ImageFileReaderType;
    typedef itk::ImageFileWriter<ImageType>         ImageFileWriterType;
    typedef itk::ImageRegionIterator<ImageType>     ImageIteratorType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>    TensorImageIOType;
    typedef itk::TensorMeshIO<ScalarType, 3, 3>     TensorMeshIOType;
    typedef TensorImageIOType::TensorImageType      TensorImageType;
    typedef itk::TensorImageToMeshFilter<TensorType,3> TensorImageToMeshType;

    bool force_reading_all = 1;
  
    bool useprolatesystem    = cl.follow(false,2,"-p","-P");
    const char*  inputfile                    = cl.follow("input.vtk",2,"-i","-I");
    const char*  referencefile                = cl.follow("referencefile.vtk",2,"-r","-r");
    const char*  prolatefile                  = cl.follow("prolate.pr",2,"-pr","-PR");
    const char*  displacementfieldfile        = cl.follow("forward.mha",2,"-f1","-F1");
    const char*  inversedisplacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
    const char*  outputcsvfile                = cl.follow("output.csv",2,"-o","-O");
    const double lambda                       = cl.follow(1.0, 2,"-l","-L");
    const double rhobegin                     = cl.follow(0.03, 2,"-rhobegin","-RHOBEGIN");
    const double rhoend                       = cl.follow(0.005, 2,"-rhoend","-RHOEND");
    const unsigned int kerneltouse            = cl.follow(0, 2,"-k","-K");
    const char* initialposition               = cl.follow("nofile", 2,"-u","-U");
  
    std::cout << "Processing optimization: " << std::endl;
    std::cout << "useprolatesystem: \t\t" << useprolatesystem << std::endl;
    std::cout << "prolatefile: \t\t\t" << prolatefile << std::endl;
    std::cout << "displacementfieldfile: \t\t" << displacementfieldfile << std::endl;
    std::cout << "inversedisplacementfieldfile: \t" << inversedisplacementfieldfile << std::endl;
    std::cout << "inputfile: \t\t\t" << inputfile << std::endl;
    std::cout << "referencefile: \t\t\t" << referencefile << std::endl;
    std::cout << "outputcsvfile: \t\t\t" <<outputcsvfile << std::endl;
    std::cout << "lambda:       \t\t\t" <<lambda << std::endl;
    std::cout << "rhobegin:     \t\t\t" <<rhobegin << std::endl;
    std::cout << "rhoend:       \t\t\t" <<rhoend << std::endl;
    std::cout << "kerneltouse:    \t\t" <<kerneltouse << std::endl;
    std::cout << "initialposition:  \t" <<initialposition << std::endl;
  
    TransformType::Pointer transform = TransformType::New();
    if (useprolatesystem || force_reading_all)
    {
      std::cout<<"reading transform : "<<prolatefile<<std::endl;
      itk::TransformFactory<TransformType>::RegisterTransform ();    
      itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
      transformreader->SetFileName( prolatefile );
      transformreader->Update();
      transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
      std::cout<<"Done."<<std::endl;
    }
  
    DisplacementFileReaderType::Pointer    displacementreader1 = DisplacementFileReaderType::New();
    DisplacementFileReaderType::Pointer    displacementreader2 = DisplacementFileReaderType::New();
    // read the displacement field images
    DisplacementFieldType::Pointer displacementfield = NULL;
    if (useprolatesystem || force_reading_all)
    {
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
    }

    DisplacementFieldType::Pointer inversedisplacementfield = NULL;
    if (useprolatesystem || force_reading_all)
    {
    
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
    }
  
    std::cout << "Reading input tensors: " << inputfile << std::endl;
    TensorMeshIOType::Pointer tensorreader1 = TensorMeshIOType::New();
    tensorreader1->SetFileName(inputfile);
    try
    {
      tensorreader1->Read();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr<<e<<std::endl;
      std::exit (EXIT_FAILURE);
    }
    std::cout<<"Done."<<std::endl;
  
    MeshType::Pointer data = tensorreader1->GetOutput();  
    
    std::cout << "Reading reference tensors: " << referencefile << std::endl;
    TensorImageIOType::Pointer tensorreader2 = TensorImageIOType::New();
    tensorreader2->SetFileName(referencefile);  
    try
    {
      tensorreader2->Read();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr<<e<<std::endl;
      std::exit (EXIT_FAILURE);
    }
    TensorImageType::Pointer reference = tensorreader2->GetOutput();
    std::cout<<"Done."<<std::endl;

    
    TensorImageType::SizeType size = reference->GetLargestPossibleRegion().GetSize();
    TensorImageType::SpacingType spacing = reference->GetSpacing();
  
    CostFunctionType::Pointer costfunction = CostFunctionType::New();
    switch(kerneltouse)
    {
	case 1:
	  costfunction->SetKernelType (CostFunctionType::BSpline);
	  break;
	case 2:
	  costfunction->SetKernelType (CostFunctionType::KaiserBessel);
	  break;
	case 0:
	default:
	  costfunction->SetKernelType (CostFunctionType::Gaussian);
	  break;
    }
  
    costfunction->SetUseProlateSpheroidalCoordinates (useprolatesystem);
    std::cout<<"setting inputs"<<std::endl;
    costfunction->SetInputs (data, reference);
  
    if (useprolatesystem || force_reading_all)
    {
      std::cout<<"setting displacement fields and transform "<<std::endl;
      costfunction->SetDisplacementFields (displacementfield, inversedisplacementfield);
      costfunction->SetTransform (transform);
    }

    CostFunctionType::Interpolator2Type::LimiterType::Pointer limiter = costfunction->GetInterpolator()->GetLimiter();
    limiter->CanineDivisionsOff();
    limiter->SetTransform (transform);
    limiter->SetDisplacementField (displacementfield);
    limiter->SetInverseDisplacementField (inversedisplacementfield);
    limiter->SetAHASegmentationType (CostFunctionType::Interpolator2Type::LimiterType::AHA_17_ZONES);
    limiter->SetVentricleSizes(15, 98 * vnl_math::pi / 180.0);
    limiter->CalculateZones();
  
    std::cout<<"defining bounds/alphas"<<std::endl;

    ParametersType alphas;
    alphas.SetSize (3 * limiter->GetNumberOfAHAZones());
    ParametersValueType bounds[2][100];
  
    OptimizerType::ScalesType scales (3 * limiter->GetNumberOfAHAZones());

    if (std::strcmp (initialposition,"nofile"))
    {
      std::cout<<"reading initial position list : "<<initialposition<<std::endl;
      std::ifstream inputliststream (initialposition);
      if(inputliststream.fail())
      {
	std::cerr << "Unable to open file: " << initialposition << std::endl;
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
	std::istringstream parse ( line );
	double* kernel = new double[3];
	for (unsigned int i=0; i<3; i++)
	  parse >> kernel[i];
	alphas[3*N+0] = kernel[0];
	alphas[3*N+1] = kernel[1];
	alphas[3*N+2] = kernel[2];
      }
    }
    
    for (unsigned int i=0; i<limiter->GetNumberOfAHAZones(); i++)
    {
      double prolatethickness = 0.233;
    
      if (useprolatesystem)
      {
	if (!std::strcmp (initialposition,"nofile"))
	{
	  alphas[3*i+0] = 2 * 0.0122985;
	  alphas[3*i+1] = 2 * 0.0736842;
	  alphas[3*i+2] = 2 * 0.118046;
	}
	
	scales[3*i+0] = 1 / prolatethickness;
	scales[3*i+1] = 1 / ( vnl_math::pi / 2.0 );
	scales[3*i+2] = 1 / ( vnl_math::pi );
      
	bounds[0][3*i+0] = 0.25 * prolatethickness / 12.0; // 0.25 mm
	bounds[0][3*i+1] = 0.5 * vnl_math::pi / 180.0;     //  0.5 degree / ~0.25 mm
	bounds[0][3*i+2] = 0.5 * vnl_math::pi / 180.0;     //  0.5 degree / ~0.25 mm
      
	bounds[1][3*i+0] = 12.0 * prolatethickness / 12.0; //  12 mm
	bounds[1][3*i+1] = vnl_math::pi / 2.0;             //  90 degree
	bounds[1][3*i+2] = vnl_math::pi;                   // 180 degree
      }
      else
      {
	
	if (!std::strcmp (initialposition,"nofile"))
	{
	  alphas[3*i+0] = 2.0 * 1.0;
	  alphas[3*i+1] = 2.0 * 1.0;
	  alphas[3*i+2] = 2.0 * 1.0;
	}
	
	scales[3*i+0] = 1/((double)(size[0]) * spacing[0]);
	scales[3*i+1] = 1/((double)(size[1]) * spacing[1]);
	scales[3*i+2] = 1/((double)(size[2]) * spacing[2]);
      
	bounds[0][3*i+0] = 0.25;      // 0.25 mm
	bounds[0][3*i+1] = 0.25;      // 0.25 mm
	bounds[0][3*i+2] = 0.25;      // 0.25 mm
      
	bounds[1][3*i+0] = 50;        // 50 mm
	bounds[1][3*i+1] = 50;        // 50 mm
	bounds[1][3*i+2] = 50;        // 50 mm
      }
    }
  
    std::cout<<"alphas : "<<std::endl<<alphas<<std::endl;
    std::cout<<"scales : "<<std::endl<<scales<<std::endl;
  
    costfunction->SetBounds (bounds);
    costfunction->SetLambda (lambda);

    costfunction->SetAttachTermType (CostFunctionType::MeanSquareDistance);
  
    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetCostFunction (costfunction);
    optimizer->SetVerbosity (0);
    optimizer->SetScales (scales);
    optimizer->SetRhoBeg (rhobegin);
    optimizer->SetRhoEnd (rhoend);
    optimizer->SetInitialPosition (alphas);
  
    try
    {
      optimizer->StartOptimization();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
    }

    std::cout<<"optimization finished"<<std::endl;

    std::cout<<"last position is "<<std::endl;

    for (unsigned int i=0; i<limiter->GetNumberOfAHAZones(); i++)
    {
      std::cout<<optimizer->GetCurrentPosition()[3*i+0]<<":"
	       <<optimizer->GetCurrentPosition()[3*i+1]<<":"
	       <<optimizer->GetCurrentPosition()[3*i+2]<<std::endl;
    }
  
    double terms[2];
    costfunction->GetOutputTerms (optimizer->GetCurrentPosition(), terms);
    std::ostringstream os;

    os << useprolatesystem << " ";
    os << prolatefile << " ";
    os << displacementfieldfile << " ";
    os << inversedisplacementfieldfile << " ";
    os << inputfile << " ";
    os << referencefile << " ";
    os << outputcsvfile << " ";
    os << lambda << " ";
    os << rhobegin << " ";
    os << rhoend << " ";
    os << " [-";
    for (unsigned int i=0; i<optimizer->GetCurrentPosition().GetSize(); i++)
    {
      os << optimizer->GetCurrentPosition()[i] <<"-";
    }
    os << terms[0] << " " << terms[1] << " ";
    os <<std::endl;

    std::cout<<"writing in output "<<outputcsvfile<<"... ";  
    std::ofstream buffer (outputcsvfile, ios::in | ios::out | ios::ate);
    buffer << os.str().c_str();
    buffer.close();
    std::cout<<"done."<<std::endl;

    return EXIT_SUCCESS;
  }

}

