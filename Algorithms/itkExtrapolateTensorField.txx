/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkExtrapolateTensorField.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_ExtrapolateTensorField_txx_
#define _itk_ExtrapolateTensorField_txx_

#include "itkExtrapolateTensorField.h"

namespace itk
{

  
  template< class TPrecision, unsigned int TDimension>
  ExtrapolateTensorField<TPrecision, TDimension>
  ::ExtrapolateTensorField()
  {
    m_Domain                   = NULL;
    m_Transform                = NULL;
    m_DisplacementField        = NULL;
    m_InverseDisplacementField = NULL;
    
    m_Warper                   = WarperType::New();
    m_DomainWarper             = WarperType::New();
    m_InverseWarper            = WarperType::New();
    m_Transformer              = TransformerType::New();
    m_DomainTransformer        = TransformerType::New();
    m_InverseTransformer       = TransformerType::New();
    m_Interpolator             = InterpolatorType::New();

    m_UseProlateCoordinates    = 0;
    
    this->SetNumberOfRequiredOutputs(1);
    this->SetNumberOfRequiredInputs(1);
  };
  
  template< class TPrecision, unsigned int TDimension>
  void
  ExtrapolateTensorField<TPrecision, TDimension>
  ::GenerateData(void)
  {
    std::cout<<"generating data"<<std::endl;
    
    if (m_UseProlateCoordinates)
    {
      m_DomainWarper->Update();
      m_DomainTransformer->Update();
      m_Warper->Update();
      m_Transformer->Update();
      m_Transformer->GetOutput()->DisconnectPipeline();
      m_DomainTransformer->GetOutput()->DisconnectPipeline();
    }
    
    m_Interpolator->Update();
    
    typename MeshType::Pointer output = m_Interpolator->GetOutput();

    if(m_UseProlateCoordinates)
    {
      output->DisconnectPipeline();
      m_InverseTransformer->SetInput (output);
      m_InverseTransformer->Update();
      m_InverseWarper->Update();
      output = m_InverseWarper->GetOutput();
    }
    
    output->DisconnectPipeline();
    
    this->GraftOutput (output);
  }

  template< class TPrecision, unsigned int TDimension>
  void
  ExtrapolateTensorField<TPrecision, TDimension>
  ::GenerateOutputInformation(void)
  {
    std::cout<<"generating information"<<std::endl;
    
    this->Superclass::GenerateOutputInformation();

    if (!this->GetInput())  itkExceptionMacro (<<"Missing input");
    if (!this->GetOutput()) itkExceptionMacro(<< "output mesh not set");
    
    typename MeshType::ConstPointer input  = this->GetInput();
    typename MeshType::Pointer      domain = this->DomainToMesh (m_Domain);

    if (m_UseProlateCoordinates)
    {
      if (!m_Transform)                itkExceptionMacro (<<"Missing transform.");
      if (!m_DisplacementField)        itkExceptionMacro (<<"Missing displacement field.");
      if (!m_InverseDisplacementField) itkExceptionMacro (<<"Missing inverse displacement field.");
      
      m_DomainWarper->SetInput (domain);
      m_DomainWarper->SetDisplacementField (m_DisplacementField);
      m_DomainWarper->SetInverseDisplacementField (m_InverseDisplacementField);
      m_DomainTransformer->SetInput (m_DomainWarper->GetOutput());
      m_DomainTransformer->SetTransform (m_Transform);
      
      m_Warper->SetInput (this->GetInput());
      m_Warper->SetDisplacementField (m_DisplacementField);
      m_Warper->SetInverseDisplacementField (m_InverseDisplacementField);
      m_Transformer->SetInput (m_Warper->GetOutput());
      m_Transformer->SetTransform (m_Transform);
      
      typename TransformType::Pointer inversetransform = TransformType::New();
      m_Transform->GetInverse (inversetransform);
      m_InverseTransformer->SetInput (m_Interpolator->GetOutput());
      m_InverseTransformer->SetTransform (inversetransform);
      m_InverseWarper->SetInput (m_InverseTransformer->GetOutput());
      m_InverseWarper->SetDisplacementField (m_InverseDisplacementField);
      m_InverseWarper->SetInverseDisplacementField (m_DisplacementField);
      
      input  = m_Transformer->GetOutput();
      domain = m_DomainTransformer->GetOutput();
    }
    
    m_Interpolator->SetInput (0, domain);
    m_Interpolator->SetInput (1, input);
    m_Interpolator->SetUsePiWorkAround(m_UseProlateCoordinates);
    
    typename MeshType::PointsContainer::Pointer    DataPoints = MeshType::PointsContainer::New();
    typename MeshType::PointDataContainer::Pointer DataPixels = MeshType::PointDataContainer::New();  

    DataPoints->Reserve (domain->GetNumberOfPoints());
    DataPixels->Reserve (domain->GetNumberOfPoints());

    this->GetOutput()->SetPoints (DataPoints);
    this->GetOutput()->SetPointData (DataPixels);
    
  }

  
  template< class TPrecision, unsigned int TDimension>
  typename ExtrapolateTensorField<TPrecision, TDimension>::MeshType::Pointer
  ExtrapolateTensorField<TPrecision, TDimension>
  ::DomainToMesh(typename ImageType::Pointer image)
  {

    typename MeshType::Pointer mesh = MeshType::New();
    
    itk::ImageRegionIterator<ImageType>   itIn(image, image->GetLargestPossibleRegion());
    typename MeshType::PointsContainer::Pointer    points = MeshType::PointsContainer::New();
    typename MeshType::PointDataContainer::Pointer data   = MeshType::PointDataContainer::New();
    mesh->SetPoints (points);
    mesh->SetPointData (data);
    
    typename MeshType::PointType x;
    unsigned int counter = 0;
    
    while(!itIn.IsAtEnd())
    {
      image->TransformIndexToPhysicalPoint (itIn.GetIndex(), x);
      if (itIn.Get() > static_cast<ScalarType>(vcl_numeric_limits<ScalarType>::epsilon()))
	mesh->SetPoint (counter++, x);
      ++itIn;
    }
    
    return mesh;
  }
  
  
  
} // end of namespace itk

#endif
