#ifndef _itk_ExtractProlateInformationCommand_h_
#define _itk_ExtractProlateInformationCommand_h_

#include "itkCommandObjectBase.h"


namespace itk {

  class ExtractProlateInformationCommand : public CommandObjectBase
  {
    
  public:
		
    typedef ExtractProlateInformationCommand Self;
    typedef CommandObjectBase Superclass;
    typedef SmartPointer <Self> Pointer;
    typedef SmartPointer <const Self> ConstPointer;
    
    itkTypeMacro(ExtractProlateInformationCommand, CommandObjectBase);
    itkNewMacro(Self);
    
    const char *GetCommandName(void)
    { return "extract"; }
    
    int Execute(int nargs, const char *args[]);
    
  protected:
    ExtractProlateInformationCommand();
    ~ExtractProlateInformationCommand();
    
  private:
    ExtractProlateInformationCommand(const Self&);
    void operator=(const Self&);
    
  };
  
}

#endif
