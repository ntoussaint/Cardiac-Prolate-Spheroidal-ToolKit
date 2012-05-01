#include "itkForwardTransformMesh2Command.h"

#include "itkProlateSpheroidalTransformTensorMeshFilter.h"

#include "itkTensorImageIO.h"
#include "itkTensorMeshIO.h"
#include <itkImageFileReader.h>
#include <itkTensorImageToMeshFilter.h>
#include <itkLimitToAHAZoneImageFilter.h>

#include "itkWarpTensorMeshFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include <itksys/SystemTools.hxx>

#include <sstream>
#include <fstream>
#include <vector>
#include <fstream>
#include <iostream>

#include "GetPot.h"
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>

namespace itk
{

  ForwardTransformMesh2Command::ForwardTransformMesh2Command()
  {
    m_ShortDescription = "Transform a tensor image or mesh in Prolate Spheroidal Coordinates";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i    [input tensor image/mesh (default : input.vtk)]\n";    
    m_LongDescription +="-pr   [prolate transform used]\n";
    m_LongDescription +="-f1   [forward displacement field (default : forward.mha)]\n";
    m_LongDescription +="-f2   [backward displacement field (default : backward.mha)]\n";
    m_LongDescription +="-o    [output vtk unstructured tensor mesh]\n";    
  }

  ForwardTransformMesh2Command::~ForwardTransformMesh2Command()
  {}

  int ForwardTransformMesh2Command::Execute (int narg, const char* arg[])
  {
    
    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout  << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char* inputfile                    = cl.follow("input.vtk",2,"-i","-I");
    const char* prolatefile                  = cl.follow("prolate.lms",2,"-pr","-PR");
    const char* outputfile                   = cl.follow("output.csv",2,"-o","-O");

    // typedefs
    typedef double                                                         ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>                           TensorImageIOType;
    typedef itk::TensorMeshIO <ScalarType, 3, 3>                           TensorMeshIOType;

    typedef TensorImageIOType::TensorImageType                             TensorImageType;
    typedef TensorImageType::PixelType                                     TensorType;

    typedef itk::Vector<ScalarType, 3>                                         DisplacementType;
    typedef itk::Image<DisplacementType, 3>                                DisplacementFieldType;
    typedef itk::ImageFileReader<DisplacementFieldType>                    DisplacementFileReaderType;

    typedef TensorMeshIOType::TensorMeshType                               MeshType;
    typedef itk::WarpTensorMeshFilter<MeshType, DisplacementFieldType>     WarperType;
    typedef itk::TensorImageToMeshFilter<TensorType, 3>                    TensorImageToMeshFilterType;

    typedef itk::ProlateSpheroidalTransformTensorMeshFilter<MeshType>      TransformerType;
    typedef TransformerType::TransformType                                 TransformType;
    typedef TransformType::InputPointType                                  PointType;

    // instantiation
    TransformType::Pointer                 transform              = TransformType::New();

    vtkDataSetReader* reader = vtkDataSetReader::New();
    reader->SetFileName (inputfile);
    reader->Update();
    vtkPointSet* pointset = vtkPointSet::SafeDownCast (reader->GetOutput());
    vtkPoints* points = pointset->GetPoints();
    vtkPoints* newpoints = vtkPoints::New();

    vtkIdType n = points->GetNumberOfPoints();
    double point[3];

    std::cout<<"reading transform "<<prolatefile<<"..."<<std::endl;
    itk::TransformFactory<TransformType>::RegisterTransform ();
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName( prolatefile );
    transformreader->Update();
  
    transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    TransformType::Pointer transform_inverse = TransformType::New();
    transform->GetInverse(transform_inverse);
    std::cout << " Done." << std::endl;

    
    for (vtkIdType i = 0; i < n; i++)
    {
      points->GetPoint(i,point);

      TransformType::PointType p; p[0] = point[0]; p[1] = point[1]; p[2] = point[2];
      
      TransformType::PointType newpoint = transform->TransformPoint (p);
      
      newpoints->InsertNextPoint(newpoint.GetDataPointer());
    }

    pointset->SetPoints (newpoints);
    vtkDataSetWriter* writer = vtkDataSetWriter::New();
    writer->SetInput (pointset);
    writer->SetFileName (outputfile);
    writer->Update();
    
    reader->Delete();
    newpoints->Delete();
    writer->Delete();
    
    return EXIT_SUCCESS;
  }
  
}
