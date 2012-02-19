#ifndef _itk_ExtractProlateInformationCommandFactory_h_
#define _itk_ExtractProlateInformationCommandFactory_h_

#include "itkObjectFactoryBase.h"

namespace itk
{
  
  class ExtractProlateInformationCommandFactory : public ObjectFactoryBase
  {
    
  public:
    typedef ExtractProlateInformationCommandFactory Self;
    typedef ObjectFactoryBase        Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;
    
    /** Class methods used to interface with the registered factories. */
    virtual const char* GetITKSourceVersion(void) const;
    virtual const char* GetDescription(void) const;
    
    /** Method for class instantiation. */
    itkFactorylessNewMacro(Self);
    static ExtractProlateInformationCommandFactory* FactoryNew() { return new ExtractProlateInformationCommandFactory;}
    
    /** Run-time type information (and related methods). */
    itkTypeMacro(ExtractProlateInformationCommandFactory, ObjectFactoryBase);
    
    /** Register one factory of this type  */
    static void RegisterOneFactory(void)
    {
      ExtractProlateInformationCommandFactory::Pointer CSFFactory = ExtractProlateInformationCommandFactory::New();
      ObjectFactoryBase::RegisterFactory( CSFFactory );
    }
    
  protected:
    ExtractProlateInformationCommandFactory();
    ~ExtractProlateInformationCommandFactory();
    
  private:
    ExtractProlateInformationCommandFactory(const Self&);
    void operator=(const Self&);
    
  };
  
}

#endif
