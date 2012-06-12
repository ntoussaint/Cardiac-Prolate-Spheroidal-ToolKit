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
#include <itkTestOptimizerCommand.h>

#include <itkNewUoaOptimizer.h>
#include <itkTestingCostFunction.h>
#include <itkPoint.h>

#include "GetPot.h"

namespace itk
{

  TestOptimizerCommand::TestOptimizerCommand()
  {
    
    m_ShortDescription = "Test the optimizer";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription += "-rhobegin [rho begin for the optimization (default is 0.03)]\n";
    m_LongDescription += "-rhoend   [rho end for the optimization (default is 0.005)]\n";
    m_LongDescription += "-x0   [starting position x]\n";
    m_LongDescription += "-y0   [starting position y]\n";
    m_LongDescription += "-z0   [starting position z]\n";
    m_LongDescription += "-sx   [scale x]\n";
    m_LongDescription += "-sy   [scale y]\n";
    m_LongDescription += "-sz   [scale z]\n";
    
  }


  TestOptimizerCommand::~TestOptimizerCommand()
  {}

  int TestOptimizerCommand::Execute (int narg, const char* arg[])
  {
  
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }

    typedef double                                  ScalarType;
    typedef Point<ScalarType, 3>                    PointType;
    
    typedef itk::NewUoaOptimizer                    OptimizerType;
    typedef OptimizerType::ParameterVectorType      ParameterVectorType;
    typedef itk::TestingCostFunction                CostFunctionType;
    typedef CostFunctionType::ParametersType        ParametersType;
    typedef CostFunctionType::ParametersValueType   ParametersValueType;
    typedef CostFunctionType::MeasureType           MeasureType;

    const double rhobegin = cl.follow(0.03,  2,"-rhobegin","-RHOBEGIN");
    const double rhoend   = cl.follow(0.005, 2,"-rhoend",  "-RHOEND");
    const double x0       = cl.follow(0,     2,"-x0",      "-X0");
    const double y0       = cl.follow(0,     2,"-y0",      "-Y0");
    const double z0       = cl.follow(0,     2,"-z0",      "-Z0");
    const double sx       = cl.follow(0,     2,"-sx",      "-SX");
    const double sy       = cl.follow(0,     2,"-sy",      "-SY");
    const double sz       = cl.follow(0,     2,"-sz",      "-SZ");
    
    std::cout << "Testing optimization: " << std::endl;
    std::cout << "rhobegin:  \t\t\t" <<rhobegin << std::endl;
    std::cout << "rhoend:    \t\t\t" <<rhoend   << std::endl;
    std::cout << "x0:      \t\t\t\t" <<x0       << std::endl;
    std::cout << "y0:      \t\t\t\t" <<y0       << std::endl;
    std::cout << "z0:      \t\t\t\t" <<z0       << std::endl;
    std::cout << "sx:      \t\t\t\t" <<sx       << std::endl;
    std::cout << "sy:      \t\t\t\t" <<sy       << std::endl;
    std::cout << "sz:      \t\t\t\t" <<sz       << std::endl;
    
    OptimizerType::Pointer optimizer = OptimizerType::New();

    ParametersType alphas; alphas.SetSize (3);
    alphas[0] = x0;
    alphas[1] = y0;
    alphas[2] = z0;    
    std::cout<<"alphas : "<<std::endl<<alphas<<std::endl;
    
    OptimizerType::ScalesType scales (3);    
    scales[0] = sx;
    scales[1] = sy;
    scales[2] = sz;    
    std::cout<<"scales : "<<std::endl<<scales<<std::endl;

    optimizer->SetVerbosity       (0);
    optimizer->SetScales          (scales);
    optimizer->SetRhoBeg          (rhobegin);
    optimizer->SetRhoEnd          (rhoend);
    optimizer->SetInitialPosition (alphas);

    CostFunctionType::Pointer costfunction = CostFunctionType::New();

    optimizer->SetCostFunction (costfunction);

    std::cout<<-vcl_numeric_limits<double>::min()<<std::endl;
    
    
    getchar();
    
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

    PointType p;
    p[0] = optimizer->GetCurrentPosition()[0];
    p[1] = optimizer->GetCurrentPosition()[1];
    p[2] = optimizer->GetCurrentPosition()[2];

    std::cout<<"last position : "<<p<<std::endl;
    
    PointType origin;
    origin[0] = 0;
    origin[1] = 0;
    origin[2] = 0;
    
    std::cout<<"\n\nlast position is distant to truth by "<< origin.EuclideanDistanceTo (p)<<" mm"<<std::endl;
  
    return EXIT_SUCCESS;
  }

}

