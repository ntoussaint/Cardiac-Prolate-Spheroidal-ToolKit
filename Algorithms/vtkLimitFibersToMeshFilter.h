/*=========================================================================

Program:   IsdImaging
Module:    $Id: vtkLimitFibersToMeshFilter.h 608 2008-01-14 08:21:23Z filus $
Language:  C++
Author:    $Author: ntoussaint $
Date:      $Date: 2008-01-14 08:21:23 +0000 (Mon, 14 Jan 2008) $
Version:   $Revision: 608 $

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _vtk_LimitFibersToMeshFilter_h_
#define _vtk_LimitFibersToMeshFilter_h_

#include "vtkPolyDataAlgorithm.h"

class vtkDelaunay3D;


class vtkLimitFibersToMeshFilter: public vtkPolyDataAlgorithm
{

 public:
  static vtkLimitFibersToMeshFilter *New();
  vtkTypeRevisionMacro(vtkLimitFibersToMeshFilter, vtkPolyDataAlgorithm);
    
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
  vtkLimitFibersToMeshFilter();
  ~vtkLimitFibersToMeshFilter();
  
  // Usual data generation method
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  

 private:
  vtkLimitFibersToMeshFilter (const vtkLimitFibersToMeshFilter&);
  void operator=(const vtkLimitFibersToMeshFilter&);

  int BooleanOperation;

  vtkDataSet* Mesh;
  vtkDelaunay3D* Mesher;
  double Bounds[6];
  
  
  
};



#endif
