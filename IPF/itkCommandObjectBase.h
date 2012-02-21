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
#ifndef _itk_CommandObjectBase_h_
#define _itk_CommandObjectBase_h_

#include "itkProcessObject.h"
/**
 */

namespace itk
{

  class CommandObjectBase : public ProcessObject
  {
  public:
    typedef CommandObjectBase        Self;
    typedef ProcessObject            Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    itkTypeMacro (CommandObjectBase, ProcessObject);

    virtual const char *GetCommandName (void) = 0;
	virtual const char *GetShortDescription (void) const;
    virtual const char *GetLongDescription (void) const;

    virtual int Execute (int nargs, const char *args[]) = 0;

  protected:
    CommandObjectBase();
    ~CommandObjectBase();
	  
	std::string m_ShortDescription;
	std::string m_LongDescription;


  private:
    CommandObjectBase (const Self&);
    void operator=(const Self&);
    
    
  };
  

} // end of namespace


#endif
