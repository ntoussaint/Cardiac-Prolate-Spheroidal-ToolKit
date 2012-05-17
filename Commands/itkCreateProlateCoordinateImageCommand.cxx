#include "itkCreateProlateCoordinateImageCommand.h"

#include "itkProlateSpheroidalTransformTensorMeshFilter.h"

#include "itkTensorImageIO.h"
#include "itkTensorMeshIO.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTensorImageToMeshFilter.h>
#include <itkLimitToAHAZoneImageFilter.h>

#include "itkWarpTensorMeshFilter.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

#include <itksys/SystemTools.hxx>

#include <vtkUnstructuredGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellType.h>

#include <sstream>
#include <fstream>
#include <vector>
#include <fstream>
#include <iostream>

#include "GetPot.h"

namespace itk
{
    
    CreateProlateCoordinateImageCommand::CreateProlateCoordinateImageCommand()
    {
        m_ShortDescription = "Simply Create an image with 3-Component Prolate Coordinates from input domain";
        m_LongDescription = m_ShortDescription;
        m_LongDescription += "\n\nUsage:\n";
        m_LongDescription +="-i    [input (tensor) image]\n";    
        m_LongDescription +="-pr   [prolate transform used]\n";
        m_LongDescription +="-f   [BACKWARD displacement field (default : backward.mha)]\n";
        m_LongDescription +="-o    [output 3-component image]\n";        
    }
    
    CreateProlateCoordinateImageCommand::~CreateProlateCoordinateImageCommand()
    {}
    
    int CreateProlateCoordinateImageCommand::Execute (int narg, const char* arg[])
    {
        
        GetPot cl(narg, const_cast<char**>(arg)); // argument parser
        if( cl.size() == 1 || cl.search(2, "--help", "-h") )
        {
            std::cout  << std::endl << this->GetLongDescription() << std::endl;
            return -1;
        }
        
        
        const char* inputfile                    = cl.follow("input.mha",2,"-i","-I");
        const char* prolatefile                  = cl.follow("prolate.tr",2,"-pr","-PR");
        const char* displacementfieldfile = cl.follow("backward.mha",2,"-f2","-F2");
        const char* outputfile                   = cl.follow("output.csv",2,"-o","-O");
        
        
        // typedefs
        typedef double                                                         ScalarType;
        typedef itk::Vector<ScalarType, 3>                                     VectorType;
        typedef itk::Image<ScalarType,3>                                       ImageType;
        typedef itk::ImageFileReader<ImageType>                                ImageFileReaderType;
        typedef itk::ImageFileWriter<ImageType>                                ImageFileWriterType;
        typedef itk::Image<VectorType, 3>                                      VectorImageType;
        typedef itk::ImageRegionIterator<ImageType>                            ImageIteratorType;
        typedef itk::ImageRegionIterator<VectorImageType>                      IteratorType;
        typedef VectorType                                                     DisplacementType;
        typedef itk::Image<DisplacementType, 3>                                DisplacementFieldType;
        typedef itk::ImageFileReader<DisplacementFieldType>                    DisplacementFileReaderType;
        typedef itk::ImageFileWriter<DisplacementFieldType>                    DisplacementFileWriterType;
        typedef itk::ProlateSpheroidalTransform<ScalarType>                    TransformType;
        typedef TransformType::InputPointType                                  PointType;
        typedef itk::LinearInterpolateImageFunction<VectorImageType, ScalarType> InterpolatorType;

        std::cout << "Reading input image: " << inputfile << std::flush;  
        ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
        reader->SetFileName(inputfile);
        try
        {
            reader->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::cout << " Done." << std::endl;
        ImageType::Pointer input = reader->GetOutput();
        
        std::cout<<"reading prolate transform"<<std::endl;
        TransformType::Pointer transform = TransformType::New();  
        itk::TransformFactory<TransformType>::RegisterTransform ();
        itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
        transformreader->SetFileName( prolatefile );
        transformreader->Update();
        transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
        std::cout << " Done." << std::endl;
        
        TransformType::Pointer inversetransform = TransformType::New();
        transform->GetInverse (inversetransform);
        
        std::cout << "Reading backward field: " << displacementfieldfile << std::flush;  
        DisplacementFileReaderType::Pointer    displacementreader = DisplacementFileReaderType::New();
        displacementreader->SetFileName(displacementfieldfile);
        try
        {
            displacementreader->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::cout << " Done." << std::endl;
        DisplacementFieldType::Pointer displacementfield = displacementreader->GetOutput();
        
        InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetInputImage (displacementfield);
        
        DisplacementFieldType::Pointer outputimage = DisplacementFieldType::New();
        outputimage->SetRegions (input->GetLargestPossibleRegion());
        outputimage->SetOrigin(input->GetOrigin());
        outputimage->SetSpacing(input->GetSpacing());  
        outputimage->SetDirection(input->GetDirection());
        outputimage->Allocate();
        outputimage->FillBuffer (0.0);
        
        IteratorType itOut(outputimage, outputimage->GetLargestPossibleRegion());
        ImageIteratorType itIn(input, input->GetLargestPossibleRegion());
        
        itk::ContinuousIndex<ScalarType, 3> index;
        PointType x;
        
        DisplacementFieldType::DirectionType inversedirection (outputimage->GetDirection().GetTranspose());
        
        while(!itOut.IsAtEnd())
        {
            outputimage->TransformIndexToPhysicalPoint (itOut.GetIndex(), x);
            bool isinside = displacementfield->TransformPhysicalPointToContinuousIndex (x, index);
            
            PointType coordinates;
            for (unsigned int i=0; i<3; i++)
                coordinates[i] = 0;
            
            if (isinside && itIn.Get())
            {
                DisplacementType d = interpolator->EvaluateAtContinuousIndex (index);
                x += d;
                
                coordinates = transform->TransformPoint(x);                
                coordinates = inversedirection * coordinates;
            }
            VectorType v;
            for (unsigned int i=0; i<3; i++)                
                v[i] = coordinates[i];
            
            itOut.Set (v);
            ++itOut; ++itIn;
        }
        
        std::cout << "Writing prolate coordinate map: " << outputfile <<"... "<<std::flush;
        DisplacementFileWriterType::Pointer writer = DisplacementFileWriterType::New();
        writer->SetFileName(outputfile);
        writer->SetInput (outputimage);
        try
        {
            writer->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::cout << " Done." << std::endl;
        
        return EXIT_SUCCESS;
    }
    
}
