#ifndef _itk_RotateProlateSpheroidCommand_h_
#define _itk_RotateProlateSpheroidCommand_h_

#include "itkCommandObjectBase.h"

namespace itk {

  class RotateProlateSpheroidCommand : public CommandObjectBase
  {
    
  public:
    
    typedef RotateProlateSpheroidCommand Self;
    typedef CommandObjectBase Superclass;
    typedef SmartPointer <Self> Pointer;
    typedef SmartPointer <const Self> ConstPointer;
    
    itkTypeMacro(RotateProlateSpheroidCommand, CommandObjectBase);
    itkNewMacro(Self);
    
    const char *GetCommandName(void)
    { return "rotate"; }
    
    int Execute(int nargs, const char *args[]);
    
  protected:
    RotateProlateSpheroidCommand();
    ~RotateProlateSpheroidCommand();
    
  private:
    RotateProlateSpheroidCommand(const Self&);
    void operator=(const Self&);
    
  };
  
}

#endif
