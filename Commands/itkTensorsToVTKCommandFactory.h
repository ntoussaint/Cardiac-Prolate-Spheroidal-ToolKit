#ifndef _itk_TensorsToVTKCommandFactory_h_
#define _itk_TensorsToVTKCommandFactory_h_

#include "itkObjectFactoryBase.h"

namespace itk
{
  
  class TensorsToVTKCommandFactory : public ObjectFactoryBase
  {
    
  public:
    typedef TensorsToVTKCommandFactory Self;
    typedef ObjectFactoryBase        Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;
    
    /** Class methods used to interface with the registered factories. */
    virtual const char* GetITKSourceVersion(void) const;
    virtual const char* GetDescription(void) const;
    
    /** Method for class instantiation. */
    itkFactorylessNewMacro(Self);
    static TensorsToVTKCommandFactory* FactoryNew() { return new TensorsToVTKCommandFactory;}
    
    /** Run-time type information (and related methods). */
    itkTypeMacro(TensorsToVTKCommandFactory, ObjectFactoryBase);
    
    /** Register one factory of this type  */
    static void RegisterOneFactory(void)
    {
      TensorsToVTKCommandFactory::Pointer CSFFactory = TensorsToVTKCommandFactory::New();
      ObjectFactoryBase::RegisterFactory( CSFFactory );
    }
    
  protected:
    TensorsToVTKCommandFactory();
    ~TensorsToVTKCommandFactory();
    
  private:
    TensorsToVTKCommandFactory(const Self&);
    void operator=(const Self&);
    
  };
  
}

#endif
