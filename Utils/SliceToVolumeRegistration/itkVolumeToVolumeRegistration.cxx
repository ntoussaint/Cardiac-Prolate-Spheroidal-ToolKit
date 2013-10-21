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
#include <itkVersorRigid3DTransform.h>

#include <itkVersorRigid3DTransformOptimizer.h>

#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkCommand.h"

#include "GetPot.h"

class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
  typedef   const OptimizerType *              OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    OptimizerPointer optimizer = 
      dynamic_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};


void PrintHelp(const char* exec)
{
  std::cout << std::endl << exec << " Usage: " << std::endl << std::endl;
  std::cout << "-f  [input fixed  3D image file]" << std::endl;
  std::cout << "-m  [input moving 3D image file]" << std::endl;
  std::cout << "-t  [optional initial 3D transform file]" << std::endl;
  std::cout << "-o  [output file base]" << std::endl;
  std::cout << "-k  [keep image container and only modify header (0/1, default: 1)]" << std::endl;
  std::cout << "-mt  [metric-type (0 : cross-correlation / 1 : mutual information / 2 : derivative-based )]" << std::endl << std::endl;
  
  std::exit(EXIT_SUCCESS);
}

int main( int argc, char *argv[] )
{
  
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
  const bool  onlychangeheader = cl.follow(1, 2,"-k","-K");
  const unsigned int metrictype = cl.follow(0, 2,"-mt","-MT");
  
  double translationthreshold = 100; 
  
  const    unsigned int    Dimension = 3;
  typedef  short   PixelType;
  
  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  
  typedef itk::VersorRigid3DTransform< double> TransformType;
  itk::TransformFactory<TransformType>::RegisterTransform ();
  
  typedef itk::VersorRigid3DTransformOptimizer                                  OptimizerType;
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
  
  MetricType::Pointer                metric = MetricType::New();
  MIMetricType::Pointer            mimetric = MIMetricType::New();
  MattesMetricType::Pointer    mattesmetric = MattesMetricType::New();
  CCMetricType::Pointer            ccmetric = CCMetricType::New();
  OptimizerType::Pointer          optimizer = OptimizerType::New();
  InterpolatorType::Pointer    interpolator = InterpolatorType::New();
  RegistrationType::Pointer    registration = RegistrationType::New();
  GaussianFilterType::Pointer      smoother = GaussianFilterType::New();
  ResampleFilterType::Pointer      resample = ResampleFilterType::New();
  TransformType::Pointer          transform = TransformType::New();

  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  
  TransformType::ParametersType startingpoint;
  startingpoint.SetSize (TransformType::ParametersDimension);
  startingpoint.Fill (0.0);

  typedef TransformType::VersorType  VersorType;
  typedef VersorType::VectorType     VectorType;
  VersorType     initialrotation;
  VectorType     axis;
  axis[0] = 0.0; axis[1] = 0.0; axis[2] = 1.0;
  const double angle = 0;
  initialrotation.Set(  axis, angle  );
  transform->SetRotation( initialrotation );
  startingpoint = transform->GetParameters();

  if (IsTransformPresent)
  {
    std::cout<<"reading initial transform : "<<transformfile<<std::endl; 
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName(transformfile);
    transformreader->Update();
    
    TransformType::Pointer startingtransform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    startingpoint = startingtransform->GetParameters();
  }
  
  unsigned int numberOfHistogramBins = 20;
  MIMetricType::HistogramType::SizeType histogramSize;
  histogramSize.SetSize(3);
  histogramSize[0] = numberOfHistogramBins;
  histogramSize[1] = numberOfHistogramBins;
  histogramSize[2] = numberOfHistogramBins;
  
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
  
  fixedImageReader->SetFileName(  fileIn1 );
  movingImageReader->SetFileName( fileIn2 );
  fixedImageReader->Update();
  movingImageReader->Update();
  FixedImageType::Pointer fixedimage = fixedImageReader->GetOutput();
  MovingImageType::Pointer movingimage = movingImageReader->GetOutput();
  
  unsigned int numberofparameters = 0;
  numberofparameters = transform->GetNumberOfParameters();
  typedef MIMetricType::ScalesType ScalesType;
  ScalesType scales( numberofparameters );
  scales.Fill( 1.0 );
  mimetric->SetDerivativeStepLengthScales(scales);
  
  FixedImageType::SizeType imagesize = fixedimage->GetLargestPossibleRegion().GetSize();
  mattesmetric->SetNumberOfSpatialSamples(imagesize[0] * imagesize[1] * imagesize[2]);
  
  registration->SetFixedImage (fixedimage);
  registration->SetMovingImage (movingimage);
  registration->SetFixedImageRegion (fixedimage->GetBufferedRegion() );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  registration->SetTransform(     transform     );
  registration->SetInitialTransformParameters( startingpoint );

  typedef OptimizerType::ScalesType OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  const double translationScale = 1.0 / 1000.0;
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  
  optimizer->SetScales( optimizerScales );
  optimizer->SetMaximumStepLength( 2.00 ); 
  optimizer->SetMinimumStepLength( 0.01 );
  optimizer->SetNumberOfIterations( 100 );
  
  switch(metrictype)
  {
      case 0:
	optimizer->MaximizeOff();
	break;
      case 1:
	optimizer->MaximizeOn();
	break;
      case 2:
	optimizer->MaximizeOn();
	break;
      default:
	optimizer->MaximizeOff();
	break;
  }  


  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  
  try 
  { 
    registration->Update();
  } 
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
  } 
  
  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
  TransformType::Pointer finaltransform = TransformType::New();
  finaltransform->SetParameters( finalParameters );
  finaltransform->SetFixedParameters( transform->GetFixedParameters() );

  TransformType::OutputVectorType offset = finaltransform->GetOffset();
  TransformType::MatrixType     rotation = finaltransform->GetMatrix();
  TransformType::Pointer inversetransform = TransformType::New();
  finaltransform->GetInverse(inversetransform);
      
  std::cout << "Resulting 3D offset = "     << offset   << std::endl;
  std::cout << "Resulting 3D rotation = \n" << rotation << std::endl;
  
  ImageType::Pointer       outputimage = ImageType::New();

  if (offset.GetNorm() > translationthreshold)
  {
    std::cerr<<"found a translation greater than "<<translationthreshold<<" mm"<<std::endl;
    std::cerr<<"The registration must have failed, skipping"<<std::endl;
    outputimage->SetRegions   (movingimage->GetLargestPossibleRegion());
    outputimage->SetOrigin    (movingimage->GetOrigin());
    outputimage->SetSpacing   (movingimage->GetSpacing());
    outputimage->SetDirection (movingimage->GetDirection());
    outputimage->Allocate();
    outputimage->SetPixelContainer (movingimage->GetPixelContainer()); 
  }
  else
  {
    if (onlychangeheader)
    {
      
      TransformType::Pointer transformtouse = inversetransform;

      typedef itk::Matrix<double,4,4> UniformMatrixType;
      UniformMatrixType matrix, transformmatrix;
      matrix.SetIdentity();
      transformmatrix.SetIdentity();

      ImageType::DirectionType direction = movingimage->GetDirection();
      ImageType::DirectionType directiontouse;
      for (unsigned int i=0; i<3; i++)
      {
        for (unsigned int j=0; j<3; j++)
	{
	  matrix[i][j] = direction[i][j];
          transformmatrix[i][j] = transformtouse->GetMatrix()[i][j];
	}
        matrix[i][3] = movingimage->GetOrigin()[i];
        transformmatrix[i][3] = transformtouse->GetOffset()[i];
      }

      UniformMatrixType matrixtouse = transformmatrix * matrix;

      ImageType::PointType origintouse;
      for (unsigned int i=0; i<3; i++)
      {
	for (unsigned int j=0; j<3; j++)
	  directiontouse[i][j] = matrixtouse[i][j];
        origintouse[i] = matrixtouse[i][3];
      }


      outputimage->SetRegions   (movingimage->GetLargestPossibleRegion());
      outputimage->SetSpacing   (movingimage->GetSpacing());
      outputimage->SetDirection (directiontouse);
      outputimage->SetOrigin    (origintouse);
      outputimage->Allocate();
      outputimage->SetPixelContainer (movingimage->GetPixelContainer()); 
    }
    else
    {
      resample->SetTransform( finaltransform );
      resample->SetInput( movingimage);
      resample->SetOutputParametersFromImage (fixedimage);
      resample->SetDefaultPixelValue( 0 );
      resample->Update();

      outputimage->SetRegions   (fixedimage->GetLargestPossibleRegion());
      outputimage->SetSpacing   (fixedimage->GetSpacing());
      outputimage->SetDirection (fixedimage->GetDirection());
      outputimage->SetOrigin    (fixedimage->GetOrigin());
      outputimage->Allocate();
      outputimage->SetPixelContainer (resample->GetOutput()->GetPixelContainer()); 
    }
  }

  std::ostringstream outputimagefile;
  outputimagefile << fileOut << ".mha";
  std::ostringstream outputtransformfile;
  outputtransformfile << fileOut << ".tfm";
  std::ostringstream outputbacktransformfile;
  outputbacktransformfile << fileOut << "-inverse.tfm";
  
  std::cout<<"writing image output image to "<<outputimagefile.str().c_str()<<std::endl;
  
  typedef itk::ImageFileWriter<ImageType> WriterFixedType;
  WriterFixedType::Pointer writer = WriterFixedType::New();
  writer->SetFileName (outputimagefile.str().c_str());
  writer->SetInput (outputimage);
  writer->Update();
  
  itk::TransformFileWriter::Pointer transformwriter = itk::TransformFileWriter::New();
  std::cout<<"writing transform to "<<outputtransformfile.str().c_str()<<std::endl;
  transformwriter->SetFileName( outputtransformfile.str().c_str() );
  transformwriter->SetInput (finaltransform);
  transformwriter->Update();

  std::cout<<"writing back transform to "<<outputbacktransformfile.str().c_str()<<std::endl;
  transformwriter->SetFileName( outputbacktransformfile.str().c_str() );
  transformwriter->SetInput (inversetransform);  
  transformwriter->Update();

  return EXIT_SUCCESS;
}
