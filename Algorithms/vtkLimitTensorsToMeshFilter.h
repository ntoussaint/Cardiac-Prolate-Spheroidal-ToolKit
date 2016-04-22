/*=========================================================================

Program:   ImagingSciences
Module:    $Id: vtkLimitTensorsToMeshFilter.h 1 2010-05-21 14:00:33Z nt08 $
Language:  C++
Author:    $Author: nt08 $
Date:      $Date: 2010-05-21 14:00:33 +0000 (Fri, 21 May 2010) $
Version:   $Revision: 1 $

Copyright (c) 2010 King's College London - Division of Imaging Sciences. All rights reserved.
See Copyright.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _vtk_LimitTensorsToMeshFilter_h_
#define _vtk_LimitTensorsToMeshFilter_h_

#include "vtkUnstructuredGridAlgorithm.h"

class vtkDelaunay3D;


class vtkLimitTensorsToMeshFilter: public vtkUnstructuredGridAlgorithm
{

 public:
  static vtkLimitTensorsToMeshFilter *New();
  vtkTypeMacro(vtkLimitTensorsToMeshFilter, vtkUnstructuredGridAlgorithm);
    
  void PrintSelf (ostream& os, vtkIndent indent){};
  
  
  /** Set the Mesh */
  void SetMesh (vtkDataSet* mesh);
  vtkGetObjectMacro (Mesh, vtkDataSet);
  

  /**
     Set the boolean operator: switch between concatenating the fibers (AND - 1)
     and removing the fibers (NOT - 0).
   */
  void SetBooleanOperation (int n)
  {
    this->BooleanOperation = n;
  }
  

  /**
     Set the boolean operator: switch between concatenating the fibers (AND - 0)
     and removing the fibers (NOT - 1).
   */
  void SetBooleanOperationToAND (void)
  {
    this->BooleanOperation = 1;
  }
  

  /**
     Set the boolean operator: switch between concatenating the fibers (AND - 0)
     and removing the fibers (NOT - 1).
   */
  void SetBooleanOperationToNOT (void)
  {
    this->BooleanOperation = 0;
  }

  vtkGetMacro (BooleanOperation, int);
  
  
 protected:
  vtkLimitTensorsToMeshFilter();
  ~vtkLimitTensorsToMeshFilter();
  
  // Usual data generation method
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  

 private:
  vtkLimitTensorsToMeshFilter (const vtkLimitTensorsToMeshFilter&);
  void operator=(const vtkLimitTensorsToMeshFilter&);

  int BooleanOperation;

  vtkDataSet* Mesh;
  vtkDelaunay3D* Mesher;
  double Bounds[6];
  
  
  
};



#endif
