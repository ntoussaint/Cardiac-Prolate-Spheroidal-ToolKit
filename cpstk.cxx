/*=========================================================================

  Program:   Prolate Spheroidal ToolKit - TTK
  Module:    $URL:  $
  Language:  C++
  Date:      $Date: 2012-02-01 16:23:47 +0000 (Wed, 09 Feb 2012) $
  Version:   $Revision: 1 $

  Copyright (c)KCL 2012. All rights reserved.
  See LICENSE.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkRotateProlateSpheroidCommandFactory.h"
#include "itkTensorsToVTKCommandFactory.h"
#include "itkTensorsToVTK2CommandFactory.h"
#include "itkCreateProlateCoordinateImageCommandFactory.h"
#include "itkExtractProlateInformationCommandFactory.h"
#include "itkTensorImageToMeshCommandFactory.h"
#include "itkTensorMeshToImageCommandFactory.h"
#include "itkExtrapolateTensorFieldCommandFactory.h"
#include "itkCreateSyntheticCardiacTensorMapCommandFactory.h"
#include "itkFindClosestProlateTransformCommandFactory.h"
#include "itkCreateProlateDomainCommandFactory.h"
#include "itkTensorMeshStatisticsCommandFactory.h"
#include "itkCumulateProlateSpheroidalDataSetsCommandFactory.h"
#include "itkLimitImageToAHAZoneCommandFactory.h"
#include "itkLimitTensorsToAHAZoneCommandFactory.h"
#include "itkApplyTransformToImageCommandFactory.h"
#include "itkApplyTransformToTensorImageCommandFactory.h"
#include "itkApplyTransformToMeshCommandFactory.h"
#include "itkResampleImage3CommandFactory.h"
#include "itkResampleTensorImage3CommandFactory.h"
#include "itkColorFibresWithProlateAngleCommandFactory.h"
#include "itkOptimizeExtrapolationKernelsCommandFactory.h"
#include "itkExtractKernelsEnveloppeCommandFactory.h"
#include "itkTensorMeshGradientCommandFactory.h"
#include "itkTensorMeshStructureCommandFactory.h"
#include "itkTensorMeshCovarianceCommandFactory.h"
#include "itkForwardTransformMeshCommandFactory.h"
#include "itkForwardTransformImageCommandFactory.h"
#include "itkBackwardTransformMeshCommandFactory.h"
#include "itkBackwardTransformImageCommandFactory.h"
#include "itkReorderDWIsCommandFactory.h"
#include "itkCropImageCommandFactory.h"
#include "itkPaintTensorImageWithAngleCommandFactory.h"

#include "itkCommandObjectFactory.h"
#include "cpstkConfigure.h"
#include "GetPot.h"

int main (int narg, char *args[])
{
  itk::RotateProlateSpheroidCommandFactory::RegisterOneFactory();
  itk::TensorsToVTKCommandFactory::RegisterOneFactory();
  itk::TensorsToVTK2CommandFactory::RegisterOneFactory();
  itk::ExtractProlateInformationCommandFactory::RegisterOneFactory();
  itk::TensorImageToMeshCommandFactory::RegisterOneFactory();
  itk::TensorMeshToImageCommandFactory::RegisterOneFactory();
  itk::ExtrapolateTensorFieldCommandFactory::RegisterOneFactory();
  itk::CreateSyntheticCardiacTensorMapCommandFactory::RegisterOneFactory();
  itk::FindClosestProlateTransformCommandFactory::RegisterOneFactory();
  itk::CreateProlateDomainCommandFactory::RegisterOneFactory();
  itk::TensorMeshStatisticsCommandFactory::RegisterOneFactory();
  itk::CumulateProlateSpheroidalDataSetsCommandFactory::RegisterOneFactory();
  itk::LimitImageToAHAZoneCommandFactory::RegisterOneFactory();
  itk::LimitTensorsToAHAZoneCommandFactory::RegisterOneFactory();
  itk::ApplyTransformToImageCommandFactory::RegisterOneFactory();
  itk::ApplyTransformToTensorImageCommandFactory::RegisterOneFactory();
  itk::ApplyTransformToMeshCommandFactory::RegisterOneFactory();
  itk::ResampleImage3CommandFactory::RegisterOneFactory();
  itk::ResampleTensorImage3CommandFactory::RegisterOneFactory();
  itk::ColorFibresWithProlateAngleCommandFactory::RegisterOneFactory();
  itk::OptimizeExtrapolationKernelsCommandFactory::RegisterOneFactory();
  itk::ExtractKernelsEnveloppeCommandFactory::RegisterOneFactory();
  itk::TensorMeshGradientCommandFactory::RegisterOneFactory();
  itk::TensorMeshStructureCommandFactory::RegisterOneFactory();
  itk::TensorMeshCovarianceCommandFactory::RegisterOneFactory();
  itk::ForwardTransformMeshCommandFactory::RegisterOneFactory();
  itk::ForwardTransformImageCommandFactory::RegisterOneFactory();
  itk::BackwardTransformMeshCommandFactory::RegisterOneFactory();
  itk::BackwardTransformImageCommandFactory::RegisterOneFactory();
  itk::CreateProlateCoordinateImageCommandFactory::RegisterOneFactory();
  itk::ReorderDWIsCommandFactory::RegisterOneFactory();
  itk::CropImageCommandFactory::RegisterOneFactory();
  itk::PaintTensorImageWithAngleCommandFactory::RegisterOneFactory();

  GetPot cl(narg, const_cast<char**>(args)); // argument parser
  if( cl.size() == 1 || cl.search(2, "--help", "-h") )
  {
    std::cout << "Software Cardiac Prolate Spheroidal ToolKit (c) 2012, version " << CPSTK_VERSION << "\n";
    std::cout << "\n";
    std::cout << "Author: Nicolas Toussaint (nicolas.toussaint@kcl.ac.uk)\n";
    std::cout << "\n";
    std::cout << "The author freely grants the non-exclusive right to use the Software for RESEARCH\n"
                 "PURPOSES ONLY. Every user of the Software will communicate to the author\n"
                 "(nicolas.toussaint@kcl.ac.uk) their remarks as to the use of the Software.\n\n";
    std::cout << "Available commands:\n";
    itk::CommandObjectFactory::PrintHelp( std::cout, 0 );
    std::cout << "\n\n\n";
    return EXIT_SUCCESS;
  }

  const char *programName = args[1];
  
  itk::CommandObjectBase::Pointer prog = itk::CommandObjectFactory::CreateCommandObject( programName );
  
  int returnValue = EXIT_FAILURE;
  
  if( !prog.IsNull() )
  {
    returnValue = prog->Execute(narg-1, (const char**)args+1);
  }
  else
  {
    std::cout << std::endl << "\'"<< programName << "\'"<< " : invalid "<< args[0] <<" command name." << std::endl;
    std::cout << "type ["<< args[0] <<" -h] for a list of commands" << std::endl << std::endl;
  }
  
  return returnValue;
  
}
