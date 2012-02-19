/*=========================================================================

  Program:   ITK Program Factory
  Module:    $RCSfile: $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) INRIA Saclay Île-de-France, Parietal Research Team. All rights reserved.
  See CodeCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCommandObjectFactory_h
#define __itkCommandObjectFactory_h

#include "itkObject.h"
#include "itkCommandObjectBase.h"

namespace itk
{

  class CommandObjectFactory : public Object
  {
  public:  
    /** Standard class typedefs. */
    typedef CommandObjectFactory           Self;
    typedef Object                   Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;
  
    /** Class Methods used to interface with the registered factories */
  
    /** Run-time type information (and related methods). */
    itkTypeMacro(CommandObjectFactory, Object);

    /** Convenient typedefs. */
    typedef ::itk::CommandObjectBase::Pointer CommandObjectBasePointer;

    /** Create the appropriate CommandObject. */
    static CommandObjectBasePointer CreateCommandObject(const char* name);

	/** Print help messages of all registered commands */
	  static void PrintHelp(std::ostream &os, Indent indent);
	  
  protected:
    CommandObjectFactory();
    ~CommandObjectFactory();
    
private:
    CommandObjectFactory(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
    
  };
  
  
} // end namespace itk

#endif
