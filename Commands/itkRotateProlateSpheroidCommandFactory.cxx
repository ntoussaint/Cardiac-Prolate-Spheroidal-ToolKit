#include "itkRotateProlateSpheroidCommandFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkRotateProlateSpheroidCommand.h"
#include "itkVersion.h"

namespace itk
{
  
  RotateProlateSpheroidCommandFactory::RotateProlateSpheroidCommandFactory()
  {
    this->RegisterOverride("itkCommandObjectBase",
			   "itkRotateProlateSpheroidCommand",
			   "Rotate a Prolate Spheroid according to a vtk file describing the antero-posterior line",
			   1,
			   CreateObjectFunction<RotateProlateSpheroidCommand>::New());
  }
  
  RotateProlateSpheroidCommandFactory::~RotateProlateSpheroidCommandFactory()
  {
  }
  
  const char* 
  RotateProlateSpheroidCommandFactory::GetITKSourceVersion(void) const
  {
    return ITK_SOURCE_VERSION;
  }
  
  const char* 
  RotateProlateSpheroidCommandFactory::GetDescription(void) const
  {
    return "Rotate a Prolate Spheroid according to a vtk file describing the antero-posterior line";
  }
  
} // end namespace itk
