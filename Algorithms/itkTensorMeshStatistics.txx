/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkTensorMeshStatistics.txx 1 2010-05-21 14:00:33Z nt08 $
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
#ifndef _itk_TensorMeshStatistics_txx_
#define _itk_TensorMeshStatistics_txx_

#include "itkTensorMeshStatistics.h"

#include <itkTensorImageIO.h>
#include <itkTensorMeshIO.h>
#include <itkCrossHelper.h>
#include <itkGradientTensorImageFilter.h>

namespace itk
{

  
  
  template< class TPrecision, unsigned int TDimension>
  void 
  TensorMeshStatistics<TPrecision, TDimension>
  ::WriteTensorMesh (typename MeshType::Pointer mesh, const char* filename) const
  {
    typename TensorMeshIO <TPrecision, TDimension>::Pointer writer = TensorMeshIO <TPrecision, TDimension>::New();
    writer->SetFileName (filename);
    writer->SetInput (mesh);
    writer->Write();
  }
  
  template< class TPrecision, unsigned int TDimension>
  void 
  TensorMeshStatistics<TPrecision, TDimension>
  ::WriteTensorImage (typename ImageType::Pointer image, const char* filename) const
  {
    typename TensorImageIO <TPrecision, TDimension>::Pointer writer = TensorImageIO <TPrecision, TDimension>::New();
    writer->SetFileName (filename);
    writer->SetInput (image);
    writer->Write();
  }
  
  template< class TPrecision, unsigned int TDimension>
  TensorMeshStatistics<TPrecision, TDimension>
  ::TensorMeshStatistics()
  {
    m_Transform                = NULL;
    m_DisplacementField        = NULL;
    m_InverseDisplacementField = NULL;
    
    m_Limiter                  = AHALimiterType::New();
    m_ImageToMeshFilter        = ImageToMeshFilterType::New();
    m_Warper                   = WarperType::New();
    m_InverseWarper            = WarperType::New();
    m_Transformer              = TransformerType::New();
    m_InverseTransformer       = TransformerType::New();

    m_StatisticsOutputType     = CovarianceMatrixNorm;
    m_CoordinatesType          = ProlateSpheroidal;
    
    m_AuxMesh                  = MeshType::New();
    m_AuxWarper                = WarperType::New();
    m_AuxTransformer           = TransformerType::New();
    m_MeshOutput               = MeshType::New();

    m_GradientFilter           = GradientFilterType::New();
    

  };
  
  template< class TPrecision, unsigned int TDimension>
  void
  TensorMeshStatistics<TPrecision, TDimension>
  ::GenerateData(void)
  {
    typename OutputImageType::Pointer output = this->GetOutput();
    typedef ImageRegionIterator<OutputImageType> OutputIteratorType;
    OutputIteratorType itOut(output, output->GetRequestedRegion());
    typename OutputImageType::IndexType idOut;
    
    typename MeshType::Pointer zone, meshtoproceed;
    typename MeshType::PointType point;
    point[0] = point[1] = point[2] = 0.0;
    unsigned int numberofzones = m_Limiter->GetNumberOfAHAZones();
    CovarianceMatrixType cov(static_cast<ScalarType>(0.0));
    TensorType tensor (static_cast<ScalarType>(0.0));
    
    for (unsigned int i=0; i<numberofzones; i++)
    {
      std::cout<<"ZONE-"<<i+1<<"::"<<std::flush;
      m_Limiter->SetAHAZone (i+1);
      std::cout<<"limiter:"<<std::flush;
      m_Limiter->Update();
      std::cout<<"imagetomesh:"<<std::flush;
      m_ImageToMeshFilter->Update();

      if (m_CoordinatesType == ProlateSpheroidal)
      {
	std::cout<<"warper:"<<std::flush;
	m_Warper->Update();
	std::cout<<"transformer:"<<std::flush;
	m_Transformer->Update();
	
	meshtoproceed = m_Transformer->GetOutput();
      }
      else
      {
	meshtoproceed = m_ImageToMeshFilter->GetOutput();
      }
      
      std::cout<<"covariance:"<<std::flush;
      cov = this->ComputeCovarianceMatrix (meshtoproceed);
      
      // std::cout<<"gradient:"<<std::flush;
      tensor = this->ComputeGradientTensor (meshtoproceed);

      zone = m_ImageToMeshFilter->GetOutput();
      // ScalarType scalar = this->EvaluateScalar (cov);
      ScalarType scalar = tensor.GetNorm();
      // ScalarType scalar = static_cast<ScalarType>(i+1);
      
      for (unsigned int j=0; j<zone->GetNumberOfPoints(); j++)
      {
	zone->GetPoint (j, &point);
	bool isinside = output->TransformPhysicalPointToIndex (point, idOut);
	if (isinside)
	{
	  itOut.SetIndex (idOut);
	  itOut.Set (scalar);
	}
      }
      
      m_MeshOutput->SetPoint (i, m_Limiter->GetZoneCentralPointCartesian());
      if (m_CoordinatesType == ProlateSpheroidal)
      {
	m_AuxMesh->SetPoint (i, m_Limiter->GetZoneCentralPointProlate());
	m_AuxMesh->SetPointData (i, tensor);
	// m_MeshOutput tensor will be calculated after the loop
      }
      else
      {
	m_MeshOutput->SetPointData (i, tensor);
      }
      
      std::cout<<"done."<<std::endl;

      std::cout<<std::endl<<"tensor :"<<std::endl;
      std::cout<<tensor<<std::endl;
   
    }
    
    if (m_CoordinatesType == ProlateSpheroidal)
    {
      m_AuxTransformer->Update();
      m_AuxWarper->Update();
    
      for (unsigned int i=0; i<numberofzones; i++)
      {
	m_AuxWarper->GetOutput()->GetPointData (i, &tensor);
	m_MeshOutput->SetPointData (i, tensor);
      }
    }
    
  }

  template< class TPrecision, unsigned int TDimension>
  typename TensorMeshStatistics<TPrecision, TDimension>::CovarianceMatrixType
  TensorMeshStatistics<TPrecision, TDimension>
  ::ComputeCovarianceMatrix (typename MeshType::Pointer mesh)
  { 
    TensorType mean (0.0);
    CovarianceMatrixType covariance3 (0.0);
    TensorType t (0.0);
    
    for (unsigned long i=0; i<mesh->GetNumberOfPoints(); i++)
    {
      mesh->GetPointData (i, &t);
      mean += t.Log();
    }
    if (mesh->GetNumberOfPoints() > 1)
      mean /= static_cast<ScalarType>(mesh->GetNumberOfPoints());
    
    CovarianceMatrixType covariance (0.0);
    
    for (unsigned long i=0; i<mesh->GetNumberOfPoints(); i++)
    {
      mesh->GetPointData (i, &t);
      covariance += CovarianceMatrixType (this->Tensor2Vec (t.Log()) - this->Tensor2Vec (mean));
    }

    if (mesh->GetNumberOfPoints() > 1)
      covariance /= static_cast<ScalarType>(mesh->GetNumberOfPoints() - 1);
    
    return covariance;
  }

  

  template< class TPrecision, unsigned int TDimension>
  typename TensorMeshStatistics<TPrecision, TDimension>::TensorType
  TensorMeshStatistics<TPrecision, TDimension>
  ::ComputeGradientTensor (typename MeshType::Pointer mesh)
  {
    typedef GradientTensorImageFilter<ImageType, ScalarType> GradientFilterType;
    typedef typename GradientFilterType::OutputImageType VectorImageType;
    
    typename MeshToImageFilterType::Pointer mesh2image = MeshToImageFilterType::New();
    mesh2image->SetInput (mesh);
    
    if (m_CoordinatesType == ProlateSpheroidal)
    {
      mesh2image->SetDomain (this->CreateDomain (mesh));
    }
    else mesh2image->SetDomain (this->CreateDomain (this->GetInput()));
    
    mesh2image->Update();

    typename ImageType::Pointer image = mesh2image->GetOutput();
    
    if (m_CoordinatesType == ProlateSpheroidal)
    {
      VectorType spacing = image->GetSpacing();
      ScalarType scalefactors[3] = {0,0,0};
      PointType p , centre; p.Fill (0.0); centre.Fill (0.0);
      for (unsigned long i=0; i<mesh->GetNumberOfPoints(); i++)
      {
	mesh->GetPoint (i, &p);
	for (unsigned int j=0; j<3; j++) centre[j] += p[j];
      }
      for (unsigned int j=0; j<3; j++) centre[j] /= static_cast<ScalarType>(mesh->GetNumberOfPoints());
      
      m_Transform->EvaluateScaleFactors (centre.GetDataPointer(), scalefactors);    
      for (unsigned int j=0; j<3; j++) spacing[j] /= scalefactors[j];
      image->DisconnectPipeline();
      image->SetSpacing (spacing);
    }

    this->WriteTensorImage (image, "tensorprolateimage.mha");
    
    typename GradientFilterType::Pointer gradientfilter = GradientFilterType::New();
    gradientfilter->SetInput (image);
    gradientfilter->Update();
    
    ImageRegionIterator<VectorImageType> it (gradientfilter->GetOutput(), gradientfilter->GetOutput()->GetLargestPossibleRegion());
    VectorType mean (0.0);
    VectorType Ox (0.0); Ox[0] = 1.0;
    unsigned long count = 0;
    while (!it.IsAtEnd())
    {
      VectorType t = it.Get();
      if (t.GetNorm())
      {
	if ( (Ox * t) < 0.0 )
	  t = -t;
	mean += t;
	count++;
      }
      ++it;
    }
    mean /= static_cast<ScalarType>(count);

    /// compute the covariance matrix of all gradients
    it.GoToBegin();
    count = 0;
    vnl_matrix_fixed<ScalarType,3,3> matrix (0.0);  
    while (!it.IsAtEnd())
    {
      VectorType t = it.Get();
      if (t.GetNorm())
      {
	if ( (Ox * t) < 0.0 )
	  t = -t;
	for (unsigned int j=0; j<3; j++)
	  for (unsigned int k=0; k<3; k++)
	    matrix[j][k] += (t[j] - mean[j])*(t[k] - mean[k]);
	count++;
      }
      ++it;
    }
    matrix /= (double)count;
    matrix *= 100;
    
    TensorType covariance; covariance.SetVnlMatrix (matrix);
    
    return covariance;
  }
  
  template< class TPrecision, unsigned int TDimension>
  typename TensorMeshStatistics<TPrecision, TDimension>::ScalarType
  TensorMeshStatistics<TPrecision, TDimension>
  ::EvaluateScalar(CovarianceMatrixType sigma)
  {
    ScalarType ret = static_cast<ScalarType>(0.0);
    switch(m_StatisticsOutputType)
    {
	case CovarianceMatrixNorm:
	  ret = static_cast<ScalarType>(std::sqrt (sigma.GetTrace()));
	  break;
	case CovarianceMatrixFA:
	  ret = static_cast<ScalarType>(sigma.GetFA());
	  break;	  
	default :
	  itkExceptionMacro (<<"StatisticsOutputType not recognized : "<<m_StatisticsOutputType);
	  break;
    }

    return ret;
  }
  
  
  template< class TPrecision, unsigned int TDimension>
  void
  TensorMeshStatistics<TPrecision, TDimension>
  ::GenerateOutputInformation(void)
  {
    this->Superclass::GenerateOutputInformation();

    if (!m_Transform)                itkExceptionMacro (<<"Missing transform.");
    if (!m_DisplacementField)        itkExceptionMacro (<<"Missing displacement field.");
    if (!m_InverseDisplacementField) itkExceptionMacro (<<"Missing inverse displacement field.");
    if (!this->GetInput())           itkExceptionMacro (<<"Missing input.");

    typename InputImageType::ConstPointer input = this->GetInput();

    m_Limiter->SetInput (input);
    m_Limiter->SetDisplacementField (m_DisplacementField);
    m_Limiter->SetInverseDisplacementField (m_InverseDisplacementField);
    m_Limiter->SetTransform (m_Transform);
    m_ImageToMeshFilter->SetInput (m_Limiter->GetOutput());
    m_Warper->SetInput (m_ImageToMeshFilter->GetOutput());
    m_Warper->SetDisplacementField (m_DisplacementField);
    m_Warper->SetInverseDisplacementField (m_InverseDisplacementField);
    m_Transformer->SetInput (m_Warper->GetOutput());
    m_Transformer->SetTransform (m_Transform);
    
    typename TransformType::Pointer inversetransform = TransformType::New();
    m_Transform->GetInverse (inversetransform);
    m_InverseTransformer->SetTransform (inversetransform);
    m_InverseWarper->SetDisplacementField (m_InverseDisplacementField);
    m_InverseWarper->SetInverseDisplacementField (m_DisplacementField);

    typename OutputImageType::Pointer output = this->GetOutput();
    output->SetRegions (input->GetLargestPossibleRegion());
    output->SetOrigin(input->GetOrigin());
    output->SetSpacing(input->GetSpacing());
    output->SetDirection(input->GetDirection());
    output->Allocate();
    output->FillBuffer (static_cast<ScalarType>(0.0));

    unsigned int numberofzones = m_Limiter->GetNumberOfAHAZones();
    
    typename MeshType::PointsContainer::Pointer    AuxDataPoints = MeshType::PointsContainer::New();
    typename MeshType::PointDataContainer::Pointer AuxDataPixels = MeshType::PointDataContainer::New();  
    AuxDataPoints->Reserve (numberofzones);
    AuxDataPixels->Reserve (numberofzones);
    m_AuxMesh->SetPoints (AuxDataPoints);
    m_AuxMesh->SetPointData (AuxDataPixels);
    m_AuxTransformer->SetInput (m_AuxMesh);
    m_AuxTransformer->SetTransform (inversetransform);
    m_AuxWarper->SetInput (m_AuxTransformer->GetOutput());
    m_AuxWarper->SetDisplacementField (m_InverseDisplacementField);
    m_AuxWarper->SetInverseDisplacementField (m_DisplacementField);
    
    typename MeshType::PointsContainer::Pointer    DataPoints = MeshType::PointsContainer::New();
    typename MeshType::PointDataContainer::Pointer DataPixels = MeshType::PointDataContainer::New();
    DataPoints->Reserve (numberofzones);
    DataPixels->Reserve (numberofzones);
    m_MeshOutput->SetPoints (DataPoints);
    m_MeshOutput->SetPointData (DataPixels);
  }
  
  template< class TPrecision, unsigned int TDimension>
  typename TensorMeshStatistics<TPrecision, TDimension>::DomainImageType::Pointer
  TensorMeshStatistics<TPrecision, TDimension>
  ::CreateDomain(typename MeshType::Pointer mesh)
  {
    double range[3][2];
    PointType p; p.Fill (0.0);
    
    for (unsigned int i=0; i<3; i++)
    {
      range[i][0] =  9999;
      range[i][1] = -9999;
    }
    
    for (unsigned long i=0; i<mesh->GetNumberOfPoints(); i++)
    {
      mesh->GetPoint (i, &p);
      
      for (unsigned int j=0; j<3; j++)
      {
	range[j][0] = std::min (p[j], range[j][0]);
	range[j][1] = std::max (p[j], range[j][1]);
      }
      
    }
    
    VectorType spacing (0.0);
    VectorType length(0.0);
    for (unsigned int i=0; i<3; i++)
    {
      length[i] = range[i][1] - range[i][0];
      spacing[i] = length[i] / std::pow ((double)(mesh->GetNumberOfPoints()), 1.0/3.0);
    }
    
    PointType origin;
    for (unsigned int k=0; k<3; k++) origin[k] = range[k][0];
    
    typename DomainImageType::SizeType size;
    typename DomainImageType::RegionType region;
    for (unsigned int i=0; i<3; i++) size[i] = (unsigned int) (length[i] / spacing[i]);
    region.SetSize (size);
    
    typename DomainImageType::Pointer image = DomainImageType::New();
    image->SetOrigin (origin);
    image->SetRegions (region);
    image->SetSpacing (spacing);
    image->Allocate();

    return image;
  }

  template< class TPrecision, unsigned int TDimension>
  typename TensorMeshStatistics<TPrecision, TDimension>::DomainImageType::Pointer
  TensorMeshStatistics<TPrecision, TDimension>
  ::CreateDomain(typename ImageType::ConstPointer input)
  {
    
    typename DomainImageType::Pointer image = DomainImageType::New();
    image->SetOrigin (input->GetOrigin());
    image->SetRegions (input->GetLargestPossibleRegion());
    image->SetSpacing (input->GetSpacing());
    image->Allocate();
    
    return image;
  }

















  template< class TPrecision, unsigned int TDimension>
  typename TensorMeshStatistics<TPrecision, TDimension>::TensorType
  TensorMeshStatistics<TPrecision, TDimension>
  ::ComputeDispersionTensor (typename MeshType::Pointer mesh)
  {
    std::cout<<std::endl<<"dispersion computation 1"<<std::endl;
    
    TensorType mean (0.0);
    TensorType t (0.0), dt (0.0);
    PointType p , centre; p.Fill (0.0); centre.Fill (0.0);
    VectorType v (0.0);
    ScalarType dp = 0.0;
    ScalarType scalefactors[3] = {0,0,0}, delta = 0.0;
    ScalarType dispersions[3] = {0,0,0};
    
    CovarianceMatrixType covariances[3];
    for (unsigned int i=0; i<3; i++) covariances[i].Fill (0.0);
    
    TensorType dispersion (0.0);
    
    // first calculate the mean tensor and central point
    for (unsigned long i=0; i<mesh->GetNumberOfPoints(); i++)
    {
      mesh->GetPointData (i, &t);
      mesh->GetPoint (i, &p);
      mean += t.Log();
      for (unsigned int j=0; j<3; j++) centre[j] += p[j];
    }
    mean /= static_cast<ScalarType>(mesh->GetNumberOfPoints());
    for (unsigned int j=0; j<3; j++) centre[j] /= static_cast<ScalarType>(mesh->GetNumberOfPoints());
    
    std::cout<<"Central Point : "<<centre<<std::endl;
    std::cout<<"Mean Tensor : "<<std::endl<<mean.Exp()<<std::endl;
    
    m_Transform->EvaluateScaleFactors (centre.GetDataPointer(), scalefactors);
    
    // iterate over the dataset
    for (unsigned long i=0; i<mesh->GetNumberOfPoints(); i++)
    {
      mesh->GetPointData (i, &t);
      mesh->GetPoint (i, &p);
      v = p - centre;
      
      // dt, a log-deviation to the mean
      dt = TensorType (t.Log() - mean);
      
      for (unsigned int j=0; j<3; j++)
      {
	
	// estimate the distance projected on direction j.
      	if (m_CoordinatesType == ProlateSpheroidal)
      	{
      	  // dpprojected = dp;
      	  // for (unsigned int k=0; k<3; k++) if (k != j) dpprojected[k] = 0.0;
      	  // distance = m_Transform->EstimateGeodesicLength2 (centre, dpprojected, 100);
      	  dp = std::abs (v[j] / scalefactors[j]);
	}
	else
	  dp = std::abs (v[j]);
	
	// |dt|^2 quantifies the log-deviation to the mean
	// |dp|^2 quantifies how far we are in the j-direction
	
	// |dt|^2 / |dp|^2 is the square of the gradient in the j direction
	// between mean and point i
	
	// delta = sqrt ( |dt|^2 / |dp|^2 ) is the gradient in the j direction
	// between mean and point i
	delta = std::sqrt ((dt * dt).GetTrace() / (dp));
	//delta = std::sqrt ( CovarianceMatrixType (this->Tensor2Vec (dt)).GetTrace() / (dp));
	
	// dispersions[j] += std::sqrt ((dt * dt).GetTrace() / dp);
	dispersions[j] += delta;
	covariances[j] += CovarianceMatrixType (this->Tensor2Vec (dt)) / (std::sqrt (dp));
	
      }
    }
    
    if (mesh->GetNumberOfPoints() > 1)
    {
      for (unsigned int i=0; i<3; i++)
      {
    	dispersions[i] /= static_cast<ScalarType>(mesh->GetNumberOfPoints() - 1);
	covariances[i] /= static_cast<ScalarType>(mesh->GetNumberOfPoints() - 1);
      }
    }    
    
    // therefore dispersions[j] is the mean gradient in the j direction
    // around the mean
    
    VectorType v1 (0.0), v2 (0.0), v3 (0.0);
    v1[0] = 1; v2[1] = 1; v3[2] = 1;
    
    dispersion =
      TensorType (dispersions[0] * v1) + 
      TensorType (dispersions[1] * v2) + 
      TensorType (dispersions[2] * v3);
  
    return (dispersion);
  }

  template< class TPrecision, unsigned int TDimension>
  typename TensorMeshStatistics<TPrecision, TDimension>::TensorType
  TensorMeshStatistics<TPrecision, TDimension>
  ::ComputeDispersionTensor2 (typename MeshType::Pointer mesh)
  {
    std::cout<<std::endl<<"dispersion computation 2"<<std::endl;
    
    TensorType mean (0.0);
    CovarianceMatrixType covariance3 (0.0);
    TensorType t1 (0.0), t2 (0.0);
    PointType p, p1, p2; p.Fill (0.0); p1.Fill (0.0); p2.Fill (0.0);
    VectorType dp (0.0);
    ScalarType distance = 0.0, squaredistance = 1.0;
    ScalarType scalefactors[3] = {0,0,0};
    
    TensorType dispersion (0.0);
    TensorType dispersions[3];
    for (unsigned int i=0; i<3; i++) dispersions[i].Fill (0.0);
    
    for (unsigned long i1=0; i1<mesh->GetNumberOfPoints(); i1++)
    {
      for (unsigned long i2=0; i2<mesh->GetNumberOfPoints(); i2++)
      {
	if (i1 == i2)
	  continue;
	
	mesh->GetPointData (i1, &t1);
	mesh->GetPointData (i2, &t2);
	mesh->GetPoint (i1, &p1);
	mesh->GetPoint (i2, &p2);
	dp = p1 - p2;
	for (unsigned j=0; j<3; j++) p[j] = (p1[j] + p2[j]) / 2.0;
	m_Transform->EvaluateScaleFactors (p.GetDataPointer(), scalefactors);

	// this represent deltaT
	dispersion = TensorType (t1.Log() - t2.Log());
	
	for (unsigned int j=0; j<3; j++)
	{
	  if (m_CoordinatesType == ProlateSpheroidal)
	  {
	    distance = scalefactors[j] / dp[j];
	  }
	  else
	    distance = std::abs (dp[j]);
	  
	  squaredistance = distance * distance;	  
	  
	  // dispersions[j] represents the gradient in direction j
	  dispersions[j] += dispersion / (2.0 * static_cast<ScalarType>(squaredistance));
	}
      }
    }

    ScalarType Nsquare = static_cast<ScalarType>(mesh->GetNumberOfPoints() - 1);
    Nsquare *= Nsquare;
    
    if (mesh->GetNumberOfPoints() > 1)
      for (unsigned int i=0; i<3; i++)
    	dispersions[i] /= Nsquare;

    VectorType v1 (0.0), v2 (0.0), v3 (0.0);
    v1[0] = 1; v2[1] = 1; v3[2] = 1;    
  
    return (dispersions[0]);
  }


        // typename TensorType::VectorType v11 (0.0), v22 (0.0),v33 (0.0), v12 (0.0), v13 (0.0),v23 (0.0);
      // if (m_CoordinatesType == ProlateSpheroidal)
      // {
      // 	v11[0] = 1.0; v22[1] = 1.0; v33[2] = 1.0;
      // }
      // else
      // {
      // 	m_Transform->EvaluateLocalBasis (m_Limiter->GetZoneCentralPointProlate(), v11, v22, v33);
      // }	

      // v23 = v22 + v33; v23 /= std::sqrt (2.0);
      // v13 = v11 + v33; v13 /= std::sqrt (2.0);
      // v12 = v11 + v22; v12 /= std::sqrt (2.0);

      // std::cout<<"v11 : "<<v11<<std::endl;
      // std::cout<<"v12 : "<<v12<<std::endl;
      // std::cout<<"v13 : "<<v13<<std::endl;
      // std::cout<<"v22 : "<<v22<<std::endl;
      // std::cout<<"v23 : "<<v23<<std::endl;
      // std::cout<<"v33 : "<<v33<<std::endl;
      
      // TensorType W1 (v11), W2 (v22), W3 (v33), W4 (v23), W5 (v13), W6 (v12);
      // W4 *= 1.0 / std::sqrt (2.0); W5 *= 1.0 / std::sqrt (2.0); W6 *= 1.0 / std::sqrt (2.0);
      
      // TensorVectorType junk;
      // junk = cov * this->Tensor2Vec (W1);
      // ScalarType eps11 = this->Tensor2Vec (W1) * junk;
      // junk = cov * this->Tensor2Vec (W2);
      // ScalarType eps22 = this->Tensor2Vec (W2) * junk;
      // junk = cov * this->Tensor2Vec (W3);
      // ScalarType eps33 = this->Tensor2Vec (W3) * junk;
      // // junk = cov * this->Tensor2Vec (W4);
      // // ScalarType eps23 = this->Tensor2Vec (W4) * junk;
      // // junk = cov * this->Tensor2Vec (W5);
      // // ScalarType eps13 = this->Tensor2Vec (W5) * junk;
      // // junk = cov * this->Tensor2Vec (W6);
      // // ScalarType eps12 = this->Tensor2Vec (W6) * junk;

      // // tensor = eps11 * W1 + eps22 * W2 + eps33 * W3 + eps23 * W4 + eps13 * W5 + eps12 * W6;
      // tensor = eps11 * W1 + eps22 * W2 + eps33 * W3;

      // tensor *= 10.0;

  
} // end of namespace itk

#endif
