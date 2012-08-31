#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkNormalVariateGenerator.h" 
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImage.h"
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkCorrelationCoefficientHistogramImageToImageMetric.h>

#include "itkTranslationTransform.h"
#include "itkAffineTransform.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkCommand.h"

#include "GetPot.h"

void PrintHelp(const char* exec)
{
  std::cout << std::endl << exec << " Usage: " << std::endl << std::endl;
  std::cout << "-f  [input fixed  2D image file]" << std::endl;
  std::cout << "-m  [input moving 2D image file]" << std::endl;
  std::cout << "-t  [optional initial 2D transform file]" << std::endl;
  std::cout << "-o  [output file base]" << std::endl;
  std::cout << "-s  [smooth the moving image (0/1, default: 1)]" << std::endl;
  std::cout << "-d  [save output transform in the 3D space (0/1, default: 0)]" << std::endl << std::endl;
  std::cout << "-m  [metric-type (0 : cross-correlation / 1 : mutual information / 2 : derivative-based )]" << std::endl << std::endl;
  
  std::exit(EXIT_SUCCESS);
}

int main( int argc, char *argv[] )
{
  
  // #undef do_affine
  // #define do_affine

  GetPot   cl(argc, argv);
  if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    PrintHelp(cl[0]);

  const bool IsInput1Present     = cl.search(2,"-f","-F");
  const bool IsInput2Present     = cl.search(2,"-m","-M");
  const bool IsTransformPresent = cl.search(2,"-t","-T");
  
  if (!IsInput1Present || !IsInput2Present)
  {
    std::cerr << "Error: Input and/or output not set." << std::endl;
    std::exit (EXIT_FAILURE);
  }
  
  const char* fileIn1        = cl.follow("NoFile",2,"-f","-f");
  const char* fileIn2        = cl.follow("NoFile",2,"-m","-M");
  const char* fileOut        = cl.follow("output",2,"-o","-O");
  const char* transformfile  = cl.follow("NoFile",2,"-t","-T");
  const bool  smoothing      = cl.follow(1, 2,"-s","-s");
  const bool  threedimension = cl.follow(0, 2,"-d","-d");
  const unsigned int metrictype = cl.follow(1, 2,"-m","-M");
  
  double translationthreshold = 20; 
  double smoothingsigma = 2;
  double initialStepLength = 4;
  
  const    unsigned int    Dimension = 2;
  const    unsigned int    Dimension3 = 3;
  typedef  short   PixelType;

  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;
  typedef itk::Image< PixelType, Dimension3 > ImageType;
  
  typedef itk::TranslationTransform< double, 2 > TransformType;
  typedef itk::TranslationTransform< double, 3 > Transform3DType;
  typedef itk::AffineTransform< double, 2 > AffineTransformType;
  itk::TransformFactory<TransformType>::RegisterTransform ();
  itk::TransformFactory<Transform3DType>::RegisterTransform ();
  itk::TransformFactory<AffineTransformType>::RegisterTransform ();

  typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
  typedef itk::MeanSquaresImageToImageMetric< FixedImageType, MovingImageType > MetricType;
  typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<FixedImageType, MovingImageType >    MIMetricType;
  typedef itk::CorrelationCoefficientHistogramImageToImageMetric<FixedImageType, MovingImageType >    CCMetricType;
  typedef itk::MattesMutualInformationImageToImageMetric<FixedImageType, MovingImageType >    MattesMetricType;
  typedef itk::LinearInterpolateImageFunction< MovingImageType, double >    InterpolatorType;
  typedef itk::ImageRegistrationMethod< FixedImageType, MovingImageType >    RegistrationType;
  typedef itk::SmoothingRecursiveGaussianImageFilter< MovingImageType, MovingImageType > GaussianFilterType;
  typedef itk::ResampleImageFilter< MovingImageType, FixedImageType > ResampleFilterType;
  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

  typedef FixedImageType::SpacingType    SpacingType;
  typedef FixedImageType::PointType      OriginType;
  typedef FixedImageType::RegionType     RegionType;
  typedef FixedImageType::SizeType       SizeType;
  typedef FixedImageType::PixelContainer PixelContainerType;

  MetricType::Pointer         metric        = MetricType::New();
  MIMetricType::Pointer     mimetric        = MIMetricType::New();
  MattesMetricType::Pointer mattesmetric    = MattesMetricType::New();
  CCMetricType::Pointer ccmetric    = CCMetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  GaussianFilterType::Pointer      smoother = GaussianFilterType::New();
  ResampleFilterType::Pointer      resample = ResampleFilterType::New();
  TransformType::Pointer              transform = TransformType::New();
  AffineTransformType::Pointer  affinetransform = AffineTransformType::New();

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  
#ifdef do_affine
  AffineTransformType::ParametersType startingpoint;
  startingpoint.SetSize (AffineTransformType::ParametersDimension);
#else
  TransformType::ParametersType startingpoint;
  startingpoint.SetSize (TransformType::ParametersDimension);
#endif 
  startingpoint.Fill (0.0);
  
  if (IsTransformPresent)
  {
    std::cout<<"reading initial transform : "<<transformfile<<std::endl; 
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName(transformfile);
    transformreader->Update();
    
#ifdef do_affine
    AffineTransformType::Pointer startingtransform = dynamic_cast<AffineTransformType*>( transformreader->GetTransformList()->front().GetPointer() );
#else
    TransformType::Pointer startingtransform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
#endif
    startingpoint = startingtransform->GetParameters();
  }
  
  unsigned int numberOfHistogramBins = 50;
  MIMetricType::HistogramType::SizeType histogramSize;
#ifdef ITK_USE_REVIEW_STATISTICS
  histogramSize.SetSize(2);
#endif
  histogramSize[0] = numberOfHistogramBins;
  histogramSize[1] = numberOfHistogramBins;
  
  mimetric->SetHistogramSize( histogramSize );
  ccmetric->SetHistogramSize( histogramSize );
  
  switch(metrictype)
  {
      case 0:
	registration->SetMetric(metric);
	break;
      case 1:
	registration->SetMetric(mimetric);
	break;
      case 2:
	registration->SetMetric(ccmetric);
	break;
      default: 
	registration->SetMetric(mimetric);
	break;
  }
  
  
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  
#ifdef do_affine
  registration->SetTransform( affinetransform );
#else
  registration->SetTransform( transform );
#endif

  fixedImageReader->SetFileName(  fileIn1 );
  movingImageReader->SetFileName( fileIn2 );
  fixedImageReader->Update();
  movingImageReader->Update();
  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  MovingImageType::Pointer movingImage = movingImageReader->GetOutput();
  FixedImageType::DirectionType iddirection;
  iddirection.SetIdentity();
  
  fixedImage->SetDirection (iddirection);
  movingImage->SetDirection (iddirection);
  
  smoother->SetInput (movingImage);
  smoother->SetSigma (smoothingsigma);
  registration->SetFixedImage (fixedImage);
  if (smoothing)
    registration->SetMovingImage (smoother->GetOutput());
  else
    registration->SetMovingImage (movingImage);
  registration->SetFixedImageRegion (fixedImage->GetBufferedRegion() );
  
  typedef OptimizerType::ScalesType OptimizerScalesType;
  unsigned int numberofparameters = 0;
  
#ifdef do_affine
  numberofparameters = affinetransform->GetNumberOfParameters();
  registration->SetInitialTransformParameters( startingpoint );
  OptimizerScalesType optimizerScales( numberofparameters );
  const double translationScale = 1.0 / 1000.0;
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = 1.0;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
#else
  numberofparameters = transform->GetNumberOfParameters();
  registration->SetInitialTransformParameters( startingpoint );
  OptimizerScalesType optimizerScales( numberofparameters );
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
#endif
  
  typedef MIMetricType::ScalesType ScalesType;
  ScalesType scales( numberofparameters );
  scales.Fill( 1.0 );
  mimetric->SetDerivativeStepLengthScales(scales);
  
  FixedImageType::SizeType imagesize = fixedImage->GetLargestPossibleRegion().GetSize();
  mattesmetric->SetNumberOfSpatialSamples(imagesize[0] * imagesize[1]);

  
  
  optimizer->SetScales( optimizerScales );
  optimizer->SetRelaxationFactor( 0.8 );
  optimizer->SetMaximumStepLength( initialStepLength ); 
  optimizer->SetMinimumStepLength( 0.01 );
  optimizer->SetNumberOfIterations( 400 );
  optimizer->MaximizeOn();
  
  try 
  { 
    registration->StartRegistration(); 
  } 
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
  } 
  
  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
  
  SpacingType offset;
  offset[0] = finalParameters[0];
  offset[1] = finalParameters[1];
  std::cout << "Resulting 2D offset = " <<offset<< std::endl;
#ifdef do_affine
  AffineTransformType::Pointer finalTransform = AffineTransformType::New();
#else
  TransformType::Pointer finalTransform = TransformType::New();
#endif
  finalTransform->SetParameters( finalParameters );
#ifdef do_affine
  finalTransform->SetFixedParameters( affinetransform->GetFixedParameters() );
#else
  finalTransform->SetFixedParameters( transform->GetFixedParameters() );
#endif  

  PixelContainerType* container;
  
  typedef itk::ImageFileReader< ImageType  > ImageReaderType;
  ImageReaderType::Pointer imagereader = ImageReaderType::New();
  imagereader->SetFileName (fileIn1);
  imagereader->Update();
  
  if (offset.GetNorm() > translationthreshold)
  {
    std::cerr<<"found a translation greater than "<<translationthreshold<<" mm"<<std::endl;
    std::cerr<<"The registration must have failed, skipping"<<std::endl; 
    container = movingImage->GetPixelContainer();
  }
  else
  {
    std::cout<<"found translation of "<<offset.GetNorm()<<" mm"<<std::endl;
    resample->SetTransform( finalTransform );
    resample->SetInput( movingImage);
    resample->SetOutputParametersFromImage (fixedImage);
    resample->SetDefaultPixelValue( 0 );
    resample->Update();
    container = resample->GetOutput()->GetPixelContainer();
  }

  std::ostringstream outputimagefile;
  outputimagefile << fileOut << ".mha";
  std::ostringstream outputtransformfile;
  outputtransformfile << fileOut << ".mat";
  
  std::cout<<"writing image output image to "<<outputimagefile.str().c_str()<<std::endl;
  ImageType::Pointer outputimage = ImageType::New();
  outputimage->SetRegions (imagereader->GetOutput()->GetLargestPossibleRegion());
  outputimage->SetOrigin (imagereader->GetOutput()->GetOrigin());
  outputimage->SetSpacing (imagereader->GetOutput()->GetSpacing());
  outputimage->SetDirection (imagereader->GetOutput()->GetDirection());
  outputimage->Allocate();
  outputimage->SetPixelContainer (container);
  typedef itk::ImageFileWriter<ImageType> WriterFixedType;
  WriterFixedType::Pointer writer = WriterFixedType::New();
  writer->SetFileName (outputimagefile.str().c_str());
  writer->SetInput (outputimage);
  writer->Update();
  
  itk::TransformFileWriter::Pointer transformwriter = itk::TransformFileWriter::New();
  std::cout<<"writing transform to "<<outputtransformfile.str().c_str()<<std::endl;
  transformwriter->SetFileName( outputtransformfile.str().c_str() );
  
  Transform3DType::Pointer transform3d = Transform3DType::New();
  ImageType::SpacingType translation;
  translation[0] = offset[0];
  translation[1] = offset[1];
  translation[2] = 0;
  ImageType::DirectionType matrix (imagereader->GetOutput()->GetDirection());    
  translation = matrix * (-translation);
  Transform3DType::ParametersType p; p.SetSize (3);
  p[0] = translation[0];
  p[1] = translation[1];
  p[2] = translation[2];
  transform3d->SetParameters (p);

  if (threedimension)
  {
    std::cout<<"3D translation : "<<translation<<std::endl;
    transformwriter->SetInput (transform3d);
  }
  else
  {
    std::cout<<"2D translation : "<<offset<<std::endl;
    transformwriter->SetInput (finalTransform);
  }
  
  transformwriter->Update();

  return EXIT_SUCCESS;
}
