/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkTensorMeshIO.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_TensorMeshIO_txx_
#define _itk_TensorMeshIO_txx_

#include "itkTensorMeshIO.h"

#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <vtksys/SystemTools.hxx>
#include <vtksys/stl/string>

namespace itk
{


  template<class T, unsigned int TensorDimension, unsigned int ImageDimension>
  bool
  TensorMeshIO<T,TensorDimension,ImageDimension>
  ::CheckExtension (const char* filename, const char* ext) const
  {
    const char* found = strstr(filename, ext);

    // return true if extension is found and is at the end of the filename
    return ( found && found == ( filename+strlen(filename)-strlen(ext) ) );
  }
  
  template<class T, unsigned int TensorDimension, unsigned int ImageDimension>
  void
  TensorMeshIO<T,TensorDimension,ImageDimension>
  ::Read (void)
  {
    
    if(m_FileName == "")
    {
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error: FileName is not set.");
    }
    
    vtksys_stl::string ext = vtksys::SystemTools::GetFilenameExtension (m_FileName.c_str());
    if( ext=="" )
    {
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error: Extension is not set.");
    }
   
    try
    {
      if ( CheckExtension(m_FileName.c_str(), ".vtk") || CheckExtension(m_FileName.c_str(), ".fib") )
      {
        this->ReadVTK (m_FileName.c_str());
      }
      else
      {
        throw itk::ExceptionObject (__FILE__,__LINE__,"Error: extension not recognized."
                                    " Supported extenstions are: .vtk and .fib");
      }
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error in TensorMeshIO::Read()");
    }
    
  }


  template<class T, unsigned int TensorDimension, unsigned int ImageDimension>
  void
  TensorMeshIO<T,TensorDimension,ImageDimension>
  ::Write (void)
  {

    if(m_Input.IsNull())
    {
      throw itk::ExceptionObject(__FILE__,__LINE__,"Error: input is not set.");
    }
    
    if(m_FileName == "")
    {
      throw itk::ExceptionObject(__FILE__,__LINE__,"Error: FileName is not set.");
    }
    
    
    vtksys_stl::string ext = vtksys::SystemTools::GetFilenameExtension (m_FileName.c_str());
    if( ext=="" )
    {
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error: Extension is not set.");
    }
    
    
    try
    {
      if ( CheckExtension(m_FileName.c_str(), ".vtk") )
      {
        this->WriteVTK (m_FileName.c_str());
      }
      else
      {
        throw itk::ExceptionObject (__FILE__,__LINE__,"Error: extension not recognized."
                                    " Supported extenstions are:\n-vtk.");
      }
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error in TensorMeshIO::Write()");
    }
    
  }

  



  template<class T, unsigned int TensorDimension, unsigned int ImageDimension>
  void
  TensorMeshIO<T,TensorDimension,ImageDimension>
  ::ReadVTK (const char* filename)
  {
    // VTK only supports 3D meshes with 3x3 tensors
    if( TensorDimension!=3 || ImageDimension!=3 )
    {
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error: VTK only supports 3D images and 3x3 tensors.");
    }
    
    typename TensorMeshType::Pointer                          output = m_Output;
    typename TensorMeshType::PointsContainer::Pointer     DataPoints = TensorMeshType::PointsContainer::New();
    typename TensorMeshType::PointDataContainer::Pointer  DataPixels = TensorMeshType::PointDataContainer::New();
    
    this->Reader->SetFileName (filename);
    this->Reader->Update();

    if (!vtkPointSet::SafeDownCast (this->Reader->GetOutput()))
    {
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error: This Reader only supports vtkPointSet subclasses.");
    }
    bool tensors_in_points = false;
    bool tensors_in_fields = false;
    
    vtkDataArray* tensors = vtkPointSet::SafeDownCast (this->Reader->GetOutput())->GetPointData()->GetTensors();
    if (!tensors) 
    {
      itkWarningMacro (<<"Warning : No Tensors in Point_Data, attempting field-data...");
      tensors = vtkPointSet::SafeDownCast (this->Reader->GetOutput())->GetPointData()->GetArray ("Tensors");
      if (!tensors)
      {
	itkWarningMacro (<<"Warning : The PointSet you are reading does not contain any tensors. \n Only point information will be passed to the output");
      }
      else
	tensors_in_fields = true;
    }
    else
      tensors_in_points = true;
    
    unsigned long NumberOfDataPoints = vtkPointSet::SafeDownCast (this->Reader->GetOutput())->GetPoints()->GetNumberOfPoints();
    DataPoints->Reserve (NumberOfDataPoints);
    DataPixels->Reserve (NumberOfDataPoints);
    
    output->SetPoints (DataPoints);
    output->SetPointData (DataPixels);
    
    for (unsigned long i=0; i<NumberOfDataPoints; i++)
    {
      double* point  = vtkPointSet::SafeDownCast (this->Reader->GetOutput())->GetPoint (i);
      typename TensorMeshType::PointType x;
      x[0] = point[0];
      x[1] = point[1];
      x[2] = point[2];
      output->SetPoint (i, x);
      
      if (tensors_in_points)
      {
	double* tensor = tensors->GetTuple (i);
	TensorType Tensor;
	Tensor[0] = tensor[0];
	Tensor[1] = tensor[1];
	Tensor[2] = tensor[4];
	Tensor[3] = tensor[2];
	Tensor[4] = tensor[5];
	Tensor[5] = tensor[8];
	output->SetPointData (i, Tensor);
      }
      else if (tensors_in_fields)
      {
	double* tensor = tensors->GetTuple (i);
	TensorType Tensor;
	Tensor[0] = tensor[0];
	Tensor[1] = tensor[1];
	Tensor[2] = tensor[2];
	Tensor[3] = tensor[3];
	Tensor[4] = tensor[4];
	Tensor[5] = tensor[5];
	output->SetPointData (i, Tensor);
      }
    }
    
  }


  
  template<class T, unsigned int TensorDimension, unsigned int ImageDimension>
  void
  TensorMeshIO<T,TensorDimension,ImageDimension>
  ::WriteVTK (const char* filename)
  {

    // VTK only supports 3D image with 3x3 tensors
    if( TensorDimension!=3 || ImageDimension!=3 )
    {
      throw itk::ExceptionObject (__FILE__,__LINE__,"Error: VTK only supports 3D images and 3x3 tensors.");
    }

    typename TensorMeshType::ConstPointer mesh = m_Input;
    
    vtkUnstructuredGrid* vtkoutput = vtkUnstructuredGrid::New();

    this->ConvertMeshToUnstructuredGrid(mesh, vtkoutput);
    
    vtkDoubleArray* data = vtkDoubleArray::New();
    data->SetNumberOfComponents (9);
    data->SetNumberOfTuples (mesh->GetNumberOfPoints());
    typename TensorMeshType::PixelType tensor (0.0);
    for (unsigned int i=0; i<mesh->GetNumberOfPoints(); i++)
    {
      mesh->GetPointData (i, &tensor);
      double vals[9];
      vals[0] = tensor[0]; vals[1] = tensor[1]; vals[2] = tensor[3];
      vals[3] = tensor[1]; vals[4] = tensor[2]; vals[5] = tensor[4];
      vals[6] = tensor[3]; vals[7] = tensor[4]; vals[8] = tensor[5];
      data->SetTuple (i, vals);
    }
    
    vtkoutput->GetPointData()->SetTensors (data);
    vtkDataSetWriter* vtkwriter = vtkDataSetWriter::New();
    vtkwriter->SetFileName (filename);
    vtkwriter->SetInput (vtkoutput);
    vtkwriter->Update();
    
    vtkoutput->Delete();
    data->Delete();
    vtkwriter->Delete();
    
  }
  
  template<class T, unsigned int TensorDimension, unsigned int ImageDimension>
  void
  TensorMeshIO<T,TensorDimension,ImageDimension>
  ::ConvertMeshToUnstructuredGrid(typename TensorMeshType::ConstPointer mesh, vtkUnstructuredGrid* vtkmesh)
  {
    typedef typename TensorMeshType::PointsContainer PointsContainer;
    
    // Get the number of points in the mesh
    int numPoints = mesh->GetNumberOfPoints();
    if(numPoints == 0)
    {
      mesh->Print(std::cerr);
      std::cerr << "no points in Grid " << std::endl;
      exit(-1);
    }
 
    // Create the vtkPoints object and set the number of points
    vtkPoints* vpoints = vtkPoints::New();
    vpoints->SetNumberOfPoints(numPoints);
    // Iterate over all the points in the itk mesh filling in
    // the vtkPoints object as we go
    typename PointsContainer::ConstPointer points = mesh->GetPoints();
    
    // In ITK the point container is not necessarily a vector, but in VTK it is
    vtkIdType VTKId = 0;
    std::map< vtkIdType, int > IndexMap;
    
    for(typename PointsContainer::ConstIterator i = points->Begin();
	i != points->End(); ++i, VTKId++)
    {
      // Get the point index from the point container iterator
      IndexMap[ VTKId ] = i->Index();
 
      // Set the vtk point at the index with the the coord array from itk
      // itk returns a const pointer, but vtk is not const correct, so
      // we have to use a const cast to get rid of the const
      vpoints->SetPoint(VTKId, const_cast<double*>(i->Value().GetDataPointer()));
    }
 
    // Set the points on the vtk grid
    vtkmesh->SetPoints(vpoints);
 
    // Setup some VTK things
    int vtkCellCount = 0; // running counter for current cell being inserted into vtk
    int numCells = mesh->GetNumberOfCells();
    int *types = new int[numCells]; // type array for vtk
    // create vtk cells and estimate the size
    vtkCellArray* cells = vtkCellArray::New();
    cells->EstimateSize(numCells, 4);

    // Setup the line visitor
    typedef CellInterfaceVisitorImplementation< TensorType, typename TensorMeshType::CellTraits, LineCell< CellInterface<typename TensorMeshType::PixelType, typename TensorMeshType::CellTraits > >, VisitVTKCellsClassType > LineVisitor;
    
    typename LineVisitor::Pointer lv =  LineVisitor::New();
    lv->SetTypeArray(types);
    lv->SetCellCounter(&vtkCellCount);
    lv->SetCellArray(cells);
    
    // Setup the triangle visitor
    typedef CellInterfaceVisitorImplementation< TensorType, typename TensorMeshType::CellTraits, TriangleCell< CellInterface<typename TensorMeshType::PixelType, typename TensorMeshType::CellTraits > >, VisitVTKCellsClassType> TriangleVisitor;
    
    typename TriangleVisitor::Pointer tv = TriangleVisitor::New();
    tv->SetTypeArray(types);
    tv->SetCellCounter(&vtkCellCount);
    tv->SetCellArray(cells);
    
    // Setup the quadrilateral visitor
    typedef CellInterfaceVisitorImplementation< TensorType, typename TensorMeshType::CellTraits, QuadrilateralCell< CellInterface<typename TensorMeshType::PixelType, typename TensorMeshType::CellTraits > >, VisitVTKCellsClassType> QuadrilateralVisitor;
    
    typename QuadrilateralVisitor::Pointer qv =  QuadrilateralVisitor::New();
    qv->SetTypeArray(types);
    qv->SetCellCounter(&vtkCellCount);
    qv->SetCellArray(cells);
    
    // Add the visitors to a multivisitor
    
    typename TensorMeshType::CellType::MultiVisitor::Pointer mv =
      TensorMeshType::CellType::MultiVisitor::New();
    
    mv->AddVisitor(tv);
    mv->AddVisitor(qv);
    mv->AddVisitor(lv);
    
    // Now ask the mesh to accept the multivisitor which
    // will Call Visit for each cell in the mesh that matches the
    // cell types of the visitors added to the MultiVisitor
    mesh->Accept(mv);
    
    // Now set the cells on the vtk grid with the type array and cell array
    vtkmesh->SetCells(types, cells);
    
    // Clean up vtk objects
    cells->Delete();
    vpoints->Delete();
    
  }
  
  
}


#endif
