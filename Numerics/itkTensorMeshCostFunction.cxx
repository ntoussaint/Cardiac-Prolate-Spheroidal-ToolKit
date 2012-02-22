/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkTensorMeshCostFunction.cxx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itkTensorMeshCostFunction_cxx
#define _itkTensorMeshCostFunction_cxx

#include "itkTensorMeshCostFunction.h"

#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkDataSetWriter.h"


namespace itk
{

TensorMeshCostFunction
::TensorMeshCostFunction()
{
  
  m_RawData          = TensorMeshType::New();
  m_InterpolatedData = TensorMeshType::New();
  m_Reference        = TensorMeshType::New();

  m_ForwardDisplacementField  = 0;
  m_BackwardDisplacementField = 0;
  m_Transform                 = 0;
  
  for (unsigned int i=0; i<100; i++)
  {
    m_Bounds[0][i] = vcl_numeric_limits<ParametersValueType>::min();
    m_Bounds[1][i] = vcl_numeric_limits<ParametersValueType>::max();
  }

  m_UseProlateSpheroidalCoordinates = 0;
  m_KernelType = TensorMeshCostFunction::Gaussian;
  
  std::cout<<"cost function OK"<<std::endl;  
}
  

TensorMeshCostFunction
::~TensorMeshCostFunction()
{
}

void   
TensorMeshCostFunction::SetInputs (TensorMeshType::Pointer data,
				   TensorMeshType::Pointer domain,
				   TensorMeshType::Pointer reference)
{
  this->Copy (data, m_RawData);
  this->Copy (domain, m_InterpolatedData);
  this->Copy (reference, m_Reference);

  this->UpdatePipeline();
}

TensorMeshCostFunction::MeasureType
TensorMeshCostFunction
::GetValue( const ParametersType & parameters ) const
{
  itkExceptionMacro (<<"SHOULD NOT USE THIS FUNCTION as it is deprecated"<<"\n"
		     <<", prefer using itkTensorMeshImageHybridCostFunction instead"<<"\n");
  
  bool Debug = false;
  
  if (!m_RawData->GetNumberOfPoints() || !m_InterpolatedData->GetNumberOfPoints() || !m_Reference->GetNumberOfPoints())
  {
    itkExceptionMacro (<<"Missing inputs : you need to feed the Raw Data, the Domain, and the Reference data through SetInputs()");
  }

  if (m_InterpolatedData->GetNumberOfPoints() != m_Reference->GetNumberOfPoints())
  {
    itkExceptionMacro (<<"The Reference Data has to contain the same points coordinates as the Domain.\n"
		       <<"whereas they have resp. "<<m_Reference->GetNumberOfPoints()<<" points and "
		       <<m_InterpolatedData->GetNumberOfPoints()<<" points.");
  }
    
  if (m_UseProlateSpheroidalCoordinates)
  {
    if ( !m_ForwardDisplacementField || !m_BackwardDisplacementField || !m_Transform )
      itkExceptionMacro (<<"Unconsistent parameters : using prolate spheroidal coord. induce the need of displacement fields and transform through SetDisplacementFields() and SetTransform()");
  }
  
  WarperType::Pointer             WarperIn;
  WarperType::Pointer             WarperOut;
  CoordinateSwitcherType::Pointer SwitcherIn;
  CoordinateSwitcherType::Pointer SwitcherOut;
  InterpolatorType::Pointer       Interpolator;
  
  SwitcherIn   = CoordinateSwitcherType::New();
  SwitcherOut  = CoordinateSwitcherType::New();
  WarperIn     = WarperType::New();
  WarperOut    = WarperType::New();
  Interpolator = InterpolatorType::New();
  
  switch(m_KernelType)
  {
      case TensorMeshCostFunction::BSpline:
	Interpolator->SetKernel (BSplineKernelFunctionType::New());
	break;
      case TensorMeshCostFunction::KaiserBessel:
	Interpolator->SetKernel (KaiserBesselKernelType::New());
	break;
      case TensorMeshCostFunction::Gaussian:
      default:
	Interpolator->SetKernel (GaussianKernelFunctionType::New());
	break;
  }
  
  std::cout<<"asking value with p = "<<parameters<<std::endl;

  std::cout<<"checking bounds of input parameters."<<std::endl;
  
  double max = 10.00;
  
  for (unsigned int i=0; i<parameters.GetSize(); i++)
  {
    if (parameters[i] < m_Bounds[0][i])
      return max;
    if (parameters[i] > m_Bounds[1][i])
      return max;
  }

  
  // copy parameters in alphas
  double alpha[3];
  for (unsigned int i=0; i<3; i++)
    alpha[i] = parameters[i];

  // instantiate input meshes
  TensorMeshType::Pointer inputdata   = TensorMeshType::New();
  TensorMeshType::Pointer inputdomain = TensorMeshType::New();

  if(Debug)
  {
    this->WriteMesh (m_RawData, "rawdata.vtk");
    this->WriteMesh (m_InterpolatedData, "interpolateddata.vtk");
    this->WriteMesh (m_Reference, "reference.vtk");
  
  }
  
  if (m_UseProlateSpheroidalCoordinates)
  {
    WarperIn->SetDisplacementField (m_ForwardDisplacementField);
    WarperIn->SetInverseDisplacementField (m_BackwardDisplacementField);
    WarperOut->SetDisplacementField (m_BackwardDisplacementField);
    WarperOut->SetInverseDisplacementField (m_ForwardDisplacementField);
    
    SwitcherIn->SetTransform (m_Transform);
    TransformType::Pointer BackTransform = TransformType::New();
    m_Transform->GetInverse (BackTransform);
    SwitcherOut->SetTransform (BackTransform);
    
    std::cout<<"warping raw data to ellipsoid"<<std::endl;
    WarperIn->SetInput (m_RawData);
    WarperIn->Update();
    std::cout<<"switching raw data to prolate"<<std::endl;
    SwitcherIn->SetInput (WarperIn->GetOutput());
    SwitcherIn->Update();
    this->Copy (SwitcherIn->GetOutput(), inputdata);
    std::cout<<"warping domain to ellipsoid"<<std::endl;
    WarperIn->SetInput (m_InterpolatedData);
    WarperIn->Update();
    std::cout<<"switching domain to prolate"<<std::endl;
    SwitcherIn->SetInput (WarperIn->GetOutput());
    SwitcherIn->Update();
    this->Copy (SwitcherIn->GetOutput(), inputdomain);
  }
  else
  {
    this->Copy (m_RawData, inputdata);
    this->Copy (m_InterpolatedData, inputdomain);
  }

  if(Debug)
  {
    this->WriteMesh (inputdata, "inputdata.vtk");
    this->WriteMesh (inputdomain, "inputdomain.vtk");
  }
  
  
  std::cout<<"interpolation..."<<std::endl;
  Interpolator->SetInput (0, inputdomain);
  Interpolator->SetInput (1, inputdata);
  Interpolator->SetUsePiWorkAround(m_UseProlateSpheroidalCoordinates);
  Interpolator->SetAlpha (alpha);  
  Interpolator->Update();
  std::cout<<"done"<<std::endl;

  TensorMeshType::Pointer interpolatoroutput = Interpolator->GetOutput();
  interpolatoroutput->DisconnectPipeline();
  
  if(Debug)
  {
    this->WriteMesh (interpolatoroutput, "interpolatoroutput.vtk");
  }
  
  TensorMeshType::Pointer outputinterpolateddata = TensorMeshType::New();
  
  if (m_UseProlateSpheroidalCoordinates)
  {
    std::cout<<"switching back interpolated data to cartesian"<<std::endl;
    SwitcherOut->SetInput (interpolatoroutput);
    SwitcherOut->Update();
    std::cout<<"warping back interpolated data out of ellipsoid"<<std::endl;
    WarperOut->SetInput (SwitcherOut->GetOutput());
    WarperOut->Update();
    this->Copy (WarperOut->GetOutput(), outputinterpolateddata);
  }
  else
  {
    this->Copy (Interpolator->GetOutput(), outputinterpolateddata);
  }

  
  if(Debug)
  {
    this->WriteMesh (outputinterpolateddata, "outputinterpolateddata.vtk");
  }
  
  
  std::cout<<"estimating distance between interpolated data and reference data"<<std::endl;
  MeasureType distance = this->EstimateMeanSquaredDistanceBetweenTensorFields (outputinterpolateddata, m_Reference);
  std::cout<<"done : "<<distance<<std::endl;
  return distance;
}


TensorMeshCostFunction::MeasureType
TensorMeshCostFunction
::EstimateMeanAngleDifferenceBetweenTensorFields( const TensorMeshType::Pointer tensors1,
						  const TensorMeshType::Pointer tensors2 ) const
{
  TensorType T1, T2, T;
  MeasureType d = 0.0, angle;
  VectorType v1, v2;
  
  PixelContainer::ConstPointer pixels1 = tensors1->GetPointData();
  PixelContainer::ConstPointer pixels2 = tensors2->GetPointData();
  PixelContainer::ConstIterator it1 = pixels1->Begin();
  PixelContainer::ConstIterator it2 = pixels2->Begin();
  
  while( (it1 != pixels1->End()) && (it2 != pixels2->End()) )
  {
    T1 = it1.Value();
    T2 = it2.Value();
    if (T1.GetTrace() > 0 && T2.GetTrace() > 0)
    {
      v1 = T1.GetEigenvector (2);
      v2 = T2.GetEigenvector (2);
      v1.Normalize();
      v2.Normalize();
      angle = std::acos (v1 * v2) * 180.0 / vnl_math::pi;      
      d += angle;
    }
  }

  d = d / (double)(tensors1->GetNumberOfPoints());
  
  return  -1.0;
}



void 
TensorMeshCostFunction
::GetDerivative( const ParametersType & parameters,
		 DerivativeType & derivative ) const
{
  /// todo !!!
} 




bool
TensorMeshCostFunction
::Copy (const TensorMeshType::Pointer input, TensorMeshType::Pointer output) const
{

  PointContainer::Pointer outputPoints = PointContainer::New();
  const PointContainer* inputPoints = input->GetPoints();

  if( inputPoints )
  {
    outputPoints->Reserve( inputPoints->Size() );
    
    PointContainer::ConstIterator inputItr = inputPoints->Begin();
    PointContainer::ConstIterator inputEnd = inputPoints->End();
    PointContainer::Iterator outputItr = outputPoints->Begin();
    
    while( inputItr != inputEnd )
    {
      outputItr.Value() = inputItr.Value();
      ++inputItr;
      ++outputItr;
    }
    
    output->SetPoints( outputPoints );
  }
  else
  {
    return false;
  }  
  
  PixelContainer::Pointer outputPointData = PixelContainer::New();
  const PixelContainer* inputPointData = input->GetPointData();

  if( inputPointData )
  {
    outputPointData->Reserve( inputPointData->Size() );

    PixelContainer::ConstIterator inputItr = inputPointData->Begin();
    PixelContainer::ConstIterator inputEnd = inputPointData->End();
    PixelContainer::Iterator outputItr = outputPointData->Begin();
    
    while( inputItr != inputEnd )
    {
      outputItr.Value() = inputItr.Value();
      ++inputItr;
      ++outputItr;
    }
    
    output->SetPointData( outputPointData );
  }
  else
  {
    return false;
  }
    
  return true;
}



bool
TensorMeshCostFunction
::WriteMesh (const TensorMeshType::Pointer mesh, const char* filename) const
{
    vtkUnstructuredGrid* output = vtkUnstructuredGrid::New();
    vtkDoubleArray*        data = vtkDoubleArray::New();
    vtkPoints*           points = vtkPoints::New();
    
    points->SetNumberOfPoints (mesh->GetNumberOfPoints());
    data->SetNumberOfComponents (9);
    data->SetNumberOfTuples (mesh->GetNumberOfPoints());

    unsigned int NN = mesh->GetNumberOfPoints();
    PointType p; p[0] = p[1] = p[2] = 0.0;
    TensorType t (0.0);
    
    for (unsigned long i=0; i<NN; i++)
    {
      mesh->GetPoint (i, &p);
      points->SetPoint (i, p.GetDataPointer());

      mesh->GetPointData (i, &t);
      double vals[9];
      vals[0] = t[0];
      vals[1] = t[1];
      vals[2] = t[3];
      vals[3] = t[1];
      vals[4] = t[2];
      vals[5] = t[4];
      vals[6] = t[3];
      vals[7] = t[4];
      vals[8] = t[5];
      
      data->SetTuple (i, vals);
    }
    
    output->SetPoints (points);
    output->GetPointData()->SetTensors (data);
    
    vtkDataSetWriter* writer = vtkDataSetWriter::New();
    writer->SetFileName (filename);
    writer->SetInput (output);
    writer->Update();

    output->Delete();
    data->Delete();
    points->Delete();
    writer->Delete();

    return true;
}





} // end namespace itk
#endif
