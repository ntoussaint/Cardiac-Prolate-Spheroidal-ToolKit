#include "itkTensorsToVTKCommandFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkTensorsToVTKCommand.h"
#include "itkVersion.h"

namespace itk
{
  
  TensorsToVTKCommandFactory::TensorsToVTKCommandFactory()
  {
    this->RegisterOverride("itkCommandObjectBase",
			   "itkTensorsToVTKCommand",
			   "Convert a tensor image (or a list of tensor images) into a vtkUnstructuredGrid structure",
			   1,
			   CreateObjectFunction<TensorsToVTKCommand>::New());
  }
  
  TensorsToVTKCommandFactory::~TensorsToVTKCommandFactory()
  {
  }
  
  const char* 
  TensorsToVTKCommandFactory::GetITKSourceVersion(void) const
  {
    return ITK_SOURCE_VERSION;
  }
  
  const char* 
  TensorsToVTKCommandFactory::GetDescription(void) const
  {
    return "Convert a tensor image (or a list of tensor images) into a vtkUnstructuredGrid structure";
  }
  
} // end namespace itk
