#include "itkRotateProlateSpheroidCommand.h"

#include <iostream>


#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"
#include "itkProlateSpheroidalTransform.h"

#include "vtkDataSetReader.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"

#include "GetPot.h"


namespace itk
{

  RotateProlateSpheroidCommand::RotateProlateSpheroidCommand()
  {
    m_ShortDescription = "Rotate a Prolate Spheroid according to a vtk file describing the antero-posterior line\n\n";
    m_LongDescription = "Usage:\n";
    m_LongDescription += "<-i input prolate file> <-a axis vtk file> <-o output prolate file>\n\n";
    m_LongDescription += m_ShortDescription;
  }

  RotateProlateSpheroidCommand::~RotateProlateSpheroidCommand()
  {}

  int RotateProlateSpheroidCommand::Execute (int narg, const char* arg[])
  {

    itk::Object::GlobalWarningDisplayOff();
    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << this->GetLongDescription() << std::endl;
      return EXIT_FAILURE;
    }

    const bool IsInputPresent = cl.search(2,"-I","-i");
    const bool IsOutputPresent = cl.search(2,"-O","-o");
    const bool IsAxisPresent = cl.search(2,"-A","-a");

  if( !IsInputPresent || !IsOutputPresent || !IsAxisPresent)
  {
    std::cerr << "Error: Input and (or) output not set." << std::endl;
    exit (-1);
  }
  
  const char* inputfile  = cl.follow("NoFile",2,"-I","-i");
  const char* axisfile   = cl.follow("NoFile",2,"-a","-A");
  const char* outputfile = cl.follow("NoFile",2,"-O","-o");
  
  typedef double ScalarType;  
  typedef itk::ProlateSpheroidalTransform<ScalarType> TransformType;  

  std::cout<<"reading transform "<<inputfile<<std::endl;
  itk::TransformFactory<TransformType>::RegisterTransform ();
  itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
  transformreader->SetFileName( inputfile );
  transformreader->Update();
  TransformType::Pointer transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
  std::cout << " Done." << std::endl;

  vtkDataSetReader* reader = vtkDataSetReader::New();
  reader->SetFileName (axisfile);
  reader->Update();
  vtkPoints* axispoints = vtkPointSet::SafeDownCast (reader->GetOutput())->GetPoints();
  
  TransformType::ParametersType parameters = transform->GetParameters();
  TransformType::PointType base, apex, septum, newseptum, c;

  unsigned int counter = 0;
  for (unsigned int i=0; i<3; i++) base[i]   = parameters[counter++];
  for (unsigned int i=0; i<3; i++) apex[i]   = parameters[counter++];
  for (unsigned int i=0; i<3; i++) septum[i] = parameters[counter++];
  TransformType::VectorType axis1, axis2, axis3, axis (0.0);
  double e2;
  
  axis1 = apex - base;
  axis2 = septum - base;
  e2 = axis2.GetNorm();
  axis1.Normalize ();
  axis2.Normalize ();
  
  axis3 = CrossProduct (axis1, axis2);
  axis2 = CrossProduct (axis3, axis1);

  for (unsigned int i=0; i<axispoints->GetNumberOfPoints(); i++)
  {
    TransformType::PointType p (axispoints->GetPoint (i));
    c = base + ( (p - base) * axis1 ) * axis1;
    TransformType::VectorType cp = p - c;
    axis += cp;
  }

  axis /= (double)(axispoints->GetNumberOfPoints());
  axis.Normalize();

  newseptum = base + e2 * axis;

  counter = 0;
  for (unsigned int i=0; i<3; i++) parameters[counter++] = base[i];
  for (unsigned int i=0; i<3; i++) parameters[counter++] = apex[i];
  for (unsigned int i=0; i<3; i++) parameters[counter++] = newseptum[i];

  transform->SetParameters (parameters);
  
  std::cout<<"writing transform "<<outputfile<<std::endl;
  itk::TransformFileWriter::Pointer transformwriter = itk::TransformFileWriter::New();
  transformwriter->SetFileName( outputfile );
  transformwriter->SetInput (transform);  
  transformwriter->Update();
  std::cout << " Done." << std::endl;

  
  return EXIT_SUCCESS;
  }
  
}
