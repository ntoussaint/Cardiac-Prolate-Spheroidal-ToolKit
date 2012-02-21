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
#include "itkCommandObjectBase.h"

namespace itk
{

  CommandObjectBase::CommandObjectBase()
  {}


  CommandObjectBase::~CommandObjectBase()
  {}


  const char *CommandObjectBase::GetShortDescription (void) const
  {
    return m_ShortDescription.c_str();
  }

  const char *CommandObjectBase::GetLongDescription (void) const
  {
    return m_LongDescription.c_str();
  }
  
}
