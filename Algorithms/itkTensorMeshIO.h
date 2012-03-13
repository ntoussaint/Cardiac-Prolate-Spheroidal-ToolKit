/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkTensorMeshIO.h 1 2010-05-21 14:00:33Z nt08 $
  Language:  C++
  Author:    $Author: nt08 $
  Date:      $Date: 2010-05-21 14:00:33 +0000 (Fri, 21 May 2010) $
  Version:   $Revision: 1 $

  Copyright (c) 2012 King's College London - Division of Imaging Sciences. All rights reserved.
  See Copyright.txt for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

  =========================================================================*/
#ifndef _itk_TensorMeshIO_h_
#define _itk_TensorMeshIO_h_

#include <itkProcessObject.h>
#include <itkMesh.h>
#include <itkTensor.h>
#include "itkLineCell.h"
#include "itkTriangleCell.h"
#include "itkQuadrilateralCell.h"
 
// VTK
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>

class vtkUnstructuredGrid;
class vtkDataSet;

namespace itk
{
  /**
     \class TensorMeshIO
     
     Author: Nicolas Toussaint
     
  */


  /**
     class to store counts of cells found in the visit pass
  */
  class CountClass
  {
  public:
    CountClass()
    {
      m_Tetra =0;
      m_QuadraticEdgeCell =0;
      m_QuadraticTriangleCellType =0;
    }
    int m_Tetra;
    int m_QuadraticEdgeCell;
    int m_QuadraticTriangleCellType;
  };
  
  template <class MeshType>
    class VisitVTKCellsClass
    {
    public:
      vtkCellArray* m_Cells;
      int* m_LastCell;
      int* m_TypeArray;
      
      // typedef the itk cells we are interested in
      typedef CellInterface< typename MeshType::PixelType, typename MeshType::CellTraits >  CellInterfaceType;
      
      typedef LineCell<CellInterfaceType> floatLineCell;
      typedef TriangleCell<CellInterfaceType>      floatTriangleCell;
      typedef QuadrilateralCell<CellInterfaceType> floatQuadrilateralCell;
      
      // Set the vtkCellArray that will be constructed
      void SetCellArray(vtkCellArray* a)
      {
	m_Cells = a;
      }
      
      // Set the cell counter pointer
      void SetCellCounter(int* i)
      {
	m_LastCell = i;
      }
      
      // Set the type array for storing the vtk cell types
      void SetTypeArray(int* i)
      {
	m_TypeArray = i;
      }
      
      // Visit a triangle and create the VTK_TRIANGLE cell
      void Visit(unsigned long, floatTriangleCell* t)
      {
	m_Cells->InsertNextCell(3,  (vtkIdType*)t->PointIdsBegin());
	m_TypeArray[*m_LastCell] = VTK_TRIANGLE;
	(*m_LastCell)++;
      }
      
      // Visit a triangle and create the VTK_QUAD cell
      void Visit(unsigned long, floatQuadrilateralCell* t)
      {
	m_Cells->InsertNextCell(4,  (vtkIdType*)t->PointIdsBegin());
	m_TypeArray[*m_LastCell] = VTK_QUAD;
	(*m_LastCell)++;
      }
      
      // Visit a line and create the VTK_LINE cell
      void Visit(unsigned long, floatLineCell* t)
      {
	m_Cells->InsertNextCell(2,  (vtkIdType*)t->PointIdsBegin());
	m_TypeArray[*m_LastCell] = VTK_LINE;
	(*m_LastCell)++;
      }
    };
  
  /**
     \Class TensorMeshIO
     This class is a ProcessObject that helps the user read and write a PointSet that contains (or not)
     some tensor information.
     It is templated with the tensor dimension and the mesh dimension. However, only implementation
     with Dimension = 3 is provided.
     The extensions allowed are .vtk and .fib.
     The Read() method will read the provided filename using a vtkDataSetReader. The dataset has to be a
     subclass of a vtkPointSet. It will pass the point and tensor information to the itk::Mesh output object.
     It is robust to the case where no tensor is present, and a warning will be thrown.

     \seealso Tensor TensorImageIO
  */
  template <class T=double, unsigned int TensorDimension=3, unsigned int MeshDimension=3>
    class ITK_EXPORT TensorMeshIO : public ProcessObject
    {
      
    public:
    
    typedef TensorMeshIO Self;
    typedef ProcessObject Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;
    
    itkNewMacro (Self);
    itkTypeMacro (TensorMeshIO, ProcessObject);
    
    /** specific typedefs */
    typedef T                                     ValueType;
    typedef Tensor<T,TensorDimension>             TensorType;
    typedef DefaultStaticMeshTraits<TensorType, TensorDimension, TensorDimension, ValueType, ValueType, TensorType> MeshTraits;
    typedef Mesh<TensorType, MeshDimension, MeshTraits> TensorMeshType;
    typedef VisitVTKCellsClass<TensorMeshType> VisitVTKCellsClassType;
    
    static const unsigned int DegreesOfFreedom = TensorDimension*(TensorDimension+1)/2;

    /** read the data */
    void Read(void);
    
    /** write the data */
    void Write(void);
    
    itkSetStringMacro(FileName);
    itkGetStringMacro(FileName);

    itkSetConstObjectMacro(Input,  TensorMeshType);
    itkGetConstObjectMacro(Input,  TensorMeshType);
    itkGetObjectMacro(Output,      TensorMeshType);

    /**
       Retrieve the vtkObject for potential user use.
       It is the output of the vtkDataSetReader.
    */
    vtkDataSet* GetVtkOutput (void)
    {
      return this->Reader->GetOutput();
    }
    
    protected:

    TensorMeshIO()
    {
      m_Output = TensorMeshType::New();
      Reader = vtkDataSetReader::New();
    }
    ~TensorMeshIO()
    {
      Reader->Delete();
    };
    
    void PrintSelf (std::ostream &os, Indent indent) const
    {
      Superclass::PrintSelf (os,indent);      
    }
  
    /** @return True if \c filename has extension \c ext */
    bool CheckExtension(const char* filename, const char* ext) const;

    /**
       This important method converts a itk::Mesh object into a vtkUnstructuredGrid.
       If the itk::Mesh has any cells they will be passed onto the vtkUnstructuredGrid as they are
       Additionnaly, the tensor information is also passed on if present.
    */
    void ConvertMeshToUnstructuredGrid(typename TensorMeshType::ConstPointer mesh, vtkUnstructuredGrid* vtkmesh);
    /**
       \todo this method is yet to be done. It could be interesting to keep track of cells
       from the vtkObject onto the itk::Mesh object.
    */
    void ConvertPointSetToMesh(vtkPointSet* vtkmesh, typename TensorMeshType::Pointer mesh);

    /**
       The only format that this class reads is vtk files. It has to be a subclass of vtkPointSet
    */
    void ReadVTK (const char* filename);
    /**
       The only format that this class reads is vtk files. It writes the object as a vtkUnstructuredGrid.
    */
    void WriteVTK (const char* filename);
    
    private:
    
    TensorMeshIO (const Self&);
    void operator=(const Self&);
    
    std::string m_FileName;
    typename TensorMeshType::ConstPointer m_Input; 
    typename TensorMeshType::Pointer m_Output;
    
    vtkDataSetReader* Reader;
    
    };
  
  

} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTensorMeshIO.txx"
#endif

#endif
