/*=========================================================================

Program:   ImagingSciences
Module:    $Id: itkSliceToVolumeRegistration.cxx 1 2010-05-21 14:00:33Z nt08 $
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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkSliceToVolumeMutualInformationCostFunction.h"
#include "itkSliceToVolumePlaneConstrainedMutualInformationCostFunction.h"
#include "itkNewUoaOptimizer.h"
#include "itkImage.h"
#include "itkTranslationTransform.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkNormalVariateGenerator.h" 
#include <itkExhaustiveOptimizer.h>
#include "itkTransformFactory.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkContinuousIndex.h"

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Parameters: " << argv[0] << std::endl;
    std::cerr << " fixedImageFile\n movingImageFile\n";
    std::cerr << " center-initialization (0/1)\n";
    std::cerr << " constrain-type (0: no constrain; 1: slice-plane constraiin; 2: respiratory-plane constrain)\n";
    std::cerr << " outputTransformfile\n [outputtranslatedimage]\n" <<std::endl;
    return EXIT_FAILURE;
    }
  
  const unsigned int Dimension3 = 3;
  const unsigned int Dimension2 = 2;
  typedef  double   PixelType;

  typedef itk::Image< PixelType, Dimension3 >  FixedImageType;
  typedef itk::Image< PixelType, Dimension2 >  MovingImageType;
  typedef itk::TranslationTransform< double, Dimension3 > TransformType;
  typedef itk::NewUoaOptimizer OptimizerType;  
  
  typedef itk::SliceToVolumeMutualInformationCostFunction<FixedImageType, MovingImageType>  NonConstrainedMetricType;
  typedef itk::SliceToVolumePlaneConstrainedMutualInformationCostFunction<FixedImageType, MovingImageType>  MetricType;
  typedef itk::ImageFileReader<FixedImageType> ImageFileReaderType;
  typedef itk::ImageFileWriter<FixedImageType> ImageFileWriterType;
  typedef itk::TransformFileWriter TransformWriterType;
  typedef itk::TranslationTransform<PixelType, Dimension3> TransformType;
  typedef TransformType::ParametersType ParametersType;
  typedef OptimizerType::ScalesType ScalesType;
  
  typedef FixedImageType::SizeType SizeType;
  typedef FixedImageType::PointType PointType;
  typedef FixedImageType::DirectionType DirectionType;
  
  typedef itk::Vector<PixelType, Dimension3> VectorType;
  
  typedef itk::ContinuousIndex<PixelType, Dimension3 > IndexType;
  
  itk::TransformFactory<TransformType>::RegisterTransform ();

  ImageFileReaderType::Pointer      reader1              = ImageFileReaderType::New();
  ImageFileReaderType::Pointer      reader2              = ImageFileReaderType::New();
  ImageFileWriterType::Pointer      writer               = ImageFileWriterType::New();
  NonConstrainedMetricType::Pointer nonconstrainedmetric = NonConstrainedMetricType::New();  
  MetricType::Pointer               metric               = MetricType::New();  
  OptimizerType::Pointer            optimizer            = OptimizerType::New();
  TransformType::Pointer            transform            = TransformType::New();
  TransformWriterType::Pointer      trwriter             = TransformWriterType::New();
  ParametersType parameters;
  ScalesType scales;
  
  double rhobeg = 5.0;
  double rhoend = 0.0000001;
  
  reader1->SetFileName (argv[1]);
  reader2->SetFileName (argv[2]);

  reader1->Update();
  reader2->Update();

  IndexType i1, i2;
  PointType p1, p2;

  FixedImageType::Pointer fixedimage = reader1->GetOutput();
  FixedImageType::Pointer movingimage = reader2->GetOutput();
  
  for (unsigned int i=0; i<Dimension3; i++)
  {
    i1[i] = (double)(fixedimage->GetLargestPossibleRegion().GetSize()[i]) / 2.0;
    i2[i] = (double)(movingimage->GetLargestPossibleRegion().GetSize()[i]) / 2.0;
  }
  
  fixedimage->TransformContinuousIndexToPhysicalPoint (i1,p1);
  movingimage->TransformContinuousIndexToPhysicalPoint (i2,p2);
  
  VectorType translation1;
  
  if (std::atoi (argv[3]) > 0)
    for (unsigned int i=0; i<Dimension3; i++)
      translation1[i] = p1[i] - p2[i];
  else
    for (unsigned int i=0; i<Dimension3; i++)
      translation1[i] = 0.0;

  std::cout<<"centering gravity points: initial translation : "<<translation1<<std::endl;  
  
  PointType origin = movingimage->GetOrigin();
  PointType neworigin = origin + translation1;

  if (std::atoi (argv[3]) > 0)
  {
    movingimage->SetOrigin (neworigin);
  }
  
  metric->SetVolumeImage (fixedimage);
  metric->SetSliceImage (movingimage);
  nonconstrainedmetric->SetVolumeImage (fixedimage);
  nonconstrainedmetric->SetSliceImage (movingimage);
  
  unsigned int constraintype = std::atoi (argv[4]);
  unsigned int paramdim = (constraintype > 0) ? Dimension2 : Dimension3;

  VectorType planeofrespiration;
  planeofrespiration[0] = 1;
  planeofrespiration[1] = 0;
  planeofrespiration[2] = 0;
  
  DirectionType respirationconstrain;
  respirationconstrain[0][0] = 0;
  respirationconstrain[1][0] = 1;
  respirationconstrain[2][0] = 0;
  respirationconstrain[0][1] = 0;
  respirationconstrain[1][1] = 0;
  respirationconstrain[2][1] = 1;
  respirationconstrain[0][2] = 1;
  respirationconstrain[1][2] = 0;
  respirationconstrain[2][2] = 0;
  
  switch (constraintype)
  {
      case 1:
	optimizer->SetCostFunction (metric);
	metric->SetDirectionOfConstrain (movingimage->GetDirection());
	break;
      case 2:
	optimizer->SetCostFunction (metric);
	metric->SetDirectionOfConstrain (respirationconstrain);
	break;
      case 0:
      default:
	optimizer->SetCostFunction (nonconstrainedmetric);
	break;
  }
  
  parameters.SetSize (paramdim);
  scales.SetSize (paramdim);
  parameters.Fill (0.0);
  scales.Fill (1.0);
  
  std::cout<<"starting position : "<<parameters<<std::endl;
  
  optimizer->SetVerbosity (0);
  optimizer->SetRhoBeg (rhobeg);
  optimizer->SetRhoEnd (rhoend);
  optimizer->SetInitialPosition (parameters);
  optimizer->SetScales (scales);
  
  optimizer->StartOptimization();
  
  for (unsigned int i=0; i<paramdim; i++)
    parameters[i] = optimizer->GetCurrentPosition()[i];
  
  std::cout<<"ending position : "<<parameters<<std::endl;

  movingimage->SetOrigin (origin);
  
  VectorType translation2;
  translation2[0] = parameters[0];
  translation2[1] = parameters[1];

  switch (constraintype)
  {
      case 1:
	translation2[2] = 0.0;
	translation2 = metric->GetDirectionOfConstrain() * translation2;
	break;
      case 2:
	translation2[2] = 0.0;
	translation2 = metric->GetDirectionOfConstrain() * translation2;
	break;
      case 0:
      default:
	translation2[2] = parameters[2];
	break;
  }
  
  ParametersType realparameters;
  realparameters.SetSize (Dimension3);
  for (unsigned int i=0; i<Dimension3; i++)
    realparameters[i] = translation1[i] + translation2[i];

  std::cout<<"extra translation found : "<<translation2<<std::endl;
  
  transform->SetParameters (realparameters);

  trwriter->SetInput (transform);
  trwriter->SetFileName (argv[5]);
  trwriter->Update();

  if (argc > 6)
  {
    FixedImageType::Pointer newimage = movingimage;
    neworigin = transform->TransformPoint (origin);
    newimage->SetOrigin (neworigin);
    
    writer->SetInput (newimage);
    writer->SetFileName (argv[6]);
    writer->Update();
  }
  
  
  return EXIT_SUCCESS;
}

