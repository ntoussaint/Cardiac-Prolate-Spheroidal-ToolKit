#ifndef _itk_RotateProlateSpheroidCommandFactory_h_
#define _itk_RotateProlateSpheroidCommandFactory_h_

#include "itkObjectFactoryBase.h"

namespace itk
{
  
  class RotateProlateSpheroidCommandFactory : public ObjectFactoryBase
  {
    
  public:
    typedef RotateProlateSpheroidCommandFactory Self;
    typedef ObjectFactoryBase        Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;
    
    /** Class methods used to interface with the registered factories. */
    virtual const char* GetITKSourceVersion(void) const;
    virtual const char* GetDescription(void) const;
    
    /** Method for class instantiation. */
    itkFactorylessNewMacro(Self);
    static RotateProlateSpheroidCommandFactory* FactoryNew() { return new RotateProlateSpheroidCommandFactory;}
    
    /** Run-time type information (and related methods). */
    itkTypeMacro(RotateProlateSpheroidCommandFactory, ObjectFactoryBase);
    
    /** Register one factory of this type  */
    static void RegisterOneFactory(void)
    {
      RotateProlateSpheroidCommandFactory::Pointer CSFFactory = RotateProlateSpheroidCommandFactory::New();
      ObjectFactoryBase::RegisterFactory( CSFFactory );
    }
    
  protected:
    RotateProlateSpheroidCommandFactory();
    ~RotateProlateSpheroidCommandFactory();
    
  private:
    RotateProlateSpheroidCommandFactory(const Self&);
    void operator=(const Self&);
    
  };
  
}

#endif
