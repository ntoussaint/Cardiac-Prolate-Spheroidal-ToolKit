#include "itkExtractProlateInformationCommandFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkExtractProlateInformationCommand.h"
#include "itkVersion.h"

namespace itk
{
  
  ExtractProlateInformationCommandFactory::ExtractProlateInformationCommandFactory()
  {
    this->RegisterOverride("itkCommandObjectBase",
			   "itkExtractProlateInformationCommand",
			   "Extract meaningful information in prolate spheroidal coordinates",
			   1,
			   CreateObjectFunction<ExtractProlateInformationCommand>::New());
  }
  
  ExtractProlateInformationCommandFactory::~ExtractProlateInformationCommandFactory()
  {
  }
  
  const char* 
  ExtractProlateInformationCommandFactory::GetITKSourceVersion(void) const
  {
    return ITK_SOURCE_VERSION;
  }
  
  const char* 
  ExtractProlateInformationCommandFactory::GetDescription(void) const
  {
    return "Extract meaningful information in prolate spheroidal coordinates";
  }
  
} // end namespace itk
