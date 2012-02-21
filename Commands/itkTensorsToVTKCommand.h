#ifndef _itk_TensorsToVTKCommand_h_
#define _itk_TensorsToVTKCommand_h_

#include "itkCommandObjectBase.h"

namespace itk {

  class TensorsToVTKCommand : public CommandObjectBase
  {
    
  public:
		
    typedef TensorsToVTKCommand Self;
    typedef CommandObjectBase Superclass;
    typedef SmartPointer <Self> Pointer;
    typedef SmartPointer <const Self> ConstPointer;
    
    itkTypeMacro(TensorsToVTKCommand, CommandObjectBase);
    itkNewMacro(Self);
    
    const char *GetCommandName(void)
    { return "itk2vtk"; }
    
    int Execute(int nargs, const char *args[]);
    
  protected:
    TensorsToVTKCommand();
    ~TensorsToVTKCommand();
    
  private:
    TensorsToVTKCommand(const Self&);
    void operator=(const Self&);
    
  };
  
}

#endif
