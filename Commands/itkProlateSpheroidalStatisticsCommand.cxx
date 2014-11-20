#include "itkProlateSpheroidalStatisticsCommand.h"

#include <itkProlateSpheroidalTransform.h>

#include <itkTransformFactory.h>
#include <itkTransformFileReader.h>
#include <itksys/SystemTools.hxx>

#include <sstream>
#include <fstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <ios>

#include "GetPot.h"


namespace itk
{

  ProlateSpheroidalStatisticsCommand::ProlateSpheroidalStatisticsCommand()
  {
    m_ShortDescription = "Extract statistics from a Prolate Spheroid";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i  [input prolate transform file]\n";    
    m_LongDescription +="-o  [output csv file containing statistics]\n";
  }

  ProlateSpheroidalStatisticsCommand::~ProlateSpheroidalStatisticsCommand()
  {}

  int ProlateSpheroidalStatisticsCommand::Execute (int narg, const char* arg[])
  {
    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout  << std::endl << this->GetLongDescription() << std::endl;
      return EXIT_FAILURE;
    }

    const bool IsInputPresent = cl.search(2,"-I","-i");
    const bool IsOutputPresent = cl.search(2,"-O","-o");

    if( !IsInputPresent || !IsOutputPresent)
    {
      std::cerr << "Error: Input and (or) output not set." << std::endl;
      exit (-1);
    }
  
    const char* inputfile  = cl.follow("NoFile",2,"-I","-i");
    const char* outputfile = cl.follow("NoFile",2,"-O","-o");
    
    typedef double ScalarType;  
    typedef ProlateSpheroidalTransform<ScalarType> TransformType;  
    typedef TransformType::PointType PointType;
    typedef TransformType::VectorType VectorType;
    
    std::cout<<"reading transform "<<inputfile<<std::endl;
    itk::TransformFactory<TransformType>::RegisterTransform ();
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( inputfile );
    transformreader->Update();
    TransformType::Pointer forward = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    std::cout << " Done." << std::endl;

    TransformType::Pointer backward = TransformType::New();
    forward->GetInverse (backward);

    TransformType::ParametersType parameters = forward->GetParameters();
    PointType base, apex, septum, midventricle_pss;
    
    unsigned int counter = 0;
    for (unsigned int i=0; i<3; i++) base[i]   = parameters[counter++];
    for (unsigned int i=0; i<3; i++) apex[i]   = parameters[counter++];
    for (unsigned int i=0; i<3; i++) septum[i] = parameters[counter++];

    midventricle_pss[0] = forward->TransformPoint (septum)[0];
    midventricle_pss[1] = vnl_math::pi / 8.0;
    midventricle_pss[2] = 0.0;    

    PointType midventricle = backward->TransformPoint (midventricle_pss);
    
    std::cout << "PSS base: \t" << forward->TransformPoint (base) << std::endl;
    std::cout << "PSS apex: \t" << forward->TransformPoint (apex) << std::endl;
    std::cout << "PSS septum: \t" << forward->TransformPoint (septum) << std::endl;
    std::cout << "PSS midventricle: \t" << forward->TransformPoint (midventricle) << std::endl;
    
    std::cout << "jac forward at septum : " <<forward->GetJacobianWithRespectToCoordinates(midventricle);
    std::cout << "jac backward at septum : " <<backward->GetJacobianWithRespectToCoordinates(midventricle_pss);
    
    
    
    float division = 0.25;
    float min = 0.25;
    float max = 10.1;
    VectorType cart, prol;
    for (unsigned int i=0; i<3; i++)
    {
      cart[i] = prol[i] = 0.0;
    }

    std::ostringstream os;
    double h[3] = {1.0, 1.0, 1.0};
    forward->EvaluateScaleFactors (midventricle_pss.GetDataPointer(), h);

    os << "mm" << " " << "xsi1" << " " << "xsi2" << " " << "xsi3" << std::endl;      
    for (float sp = min; sp < max; sp += division)
    {
      cart[0] = cart[1] = cart[2] = sp;
      std::cout << " cart vector: " << cart << std::endl;  
      for (unsigned int i = 0; i<3; i++)
	prol[i] = h[i] * cart[i];
      std::cout << " prol vector: " << prol << std::endl;  
      os << sp << " " << prol[0] << " " << prol[1] << " " << prol[2] << std::endl;      
    }
    
    std::cout << "Writing output csv file: " << outputfile << " ... " << std::flush;

    std::ofstream buffer (outputfile, std::ios::out);
    buffer << os.str().c_str();
    buffer.close();  

    std::cout<<"done."<<std::endl;
    
    return EXIT_SUCCESS;
  }
  
}
