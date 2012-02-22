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
#include "itkTensorsToVTKCommandFactory.h"
#include "itkRotateProlateSpheroidCommandFactory.h"
#include "itkExtractProlateInformationCommandFactory.h"
#include "itkTensorImageToMeshCommandFactory.h"
#include "itkTensorMeshToImageCommandFactory.h"
#include "itkExtrapolateTensorFieldCommandFactory.h"
#include "itkCreateSyntheticCardiacTensorMapCommandFactory.h"
#include "itkFindClosestProlateTransformCommandFactory.h"
#include "itkCreateProlateDomainCommandFactory.h"
#include "itkTensorMeshStatisticsCommandFactory.h"
#include "itkCreateProlateAtlasCommandFactory.h"
#include "itkLimitImageToAHAZoneCommandFactory.h"
#include "itkLimitTensorImageToAHAZoneCommandFactory.h"

#include "itkCommandObjectFactory.h"

#include "ptkConfigure.h"

int main (int narg, char *args[])
{

  itk::TensorsToVTKCommandFactory::RegisterOneFactory();
  itk::RotateProlateSpheroidCommandFactory::RegisterOneFactory();
  itk::ExtractProlateInformationCommandFactory::RegisterOneFactory();
  itk::TensorImageToMeshCommandFactory::RegisterOneFactory();
  itk::TensorMeshToImageCommandFactory::RegisterOneFactory();
  itk::ExtrapolateTensorFieldCommandFactory::RegisterOneFactory();
  itk::CreateSyntheticCardiacTensorMapCommandFactory::RegisterOneFactory();
  itk::FindClosestProlateTransformCommandFactory::RegisterOneFactory();
  itk::CreateProlateDomainCommandFactory::RegisterOneFactory();
  itk::TensorMeshStatisticsCommandFactory::RegisterOneFactory();
  itk::CreateProlateAtlasCommandFactory::RegisterOneFactory();
  itk::LimitImageToAHAZoneCommandFactory::RegisterOneFactory();
  itk::LimitTensorImageToAHAZoneCommandFactory::RegisterOneFactory();
  
  if (narg<2) {
    std::cout << "Software Prolate Spheroidal ToolKit (c)KCL 2012, version " << PTK_VERSION << "\n";
    std::cout << "\n";
    std::cout << "Author: Nicolas Toussaint (nicolas.toussaint@kcl.ac.uk)\n";
    std::cout << "\n";
    std::cout << "King's College London freely grants the non-exclusive right to use the Software for RESEARCH\n"
                 "PURPOSES ONLY. Every user of the Software will communicate to the author\n"
                 "(nicolas.toussaint@kcl.ac.uk) their remarks as to the use of the Software.\n\n";
    std::cout << "Available commands:\n";
    itk::CommandObjectFactory::PrintHelp( std::cout, 0 );
    std::cout << "\n\n\n";
    return EXIT_FAILURE;
  }

  const char *programName = args[1];
  
  itk::CommandObjectBase::Pointer prog = itk::CommandObjectFactory::CreateCommandObject( programName );
  
  int returnValue = EXIT_SUCCESS;
  
  if( !prog.IsNull() )
  {
    returnValue = prog->Execute(narg-1, (const char**)args+1);
  }
  else {
    std::cout << "Command not found" << std::endl;
  }
  
  return returnValue;
  
}
