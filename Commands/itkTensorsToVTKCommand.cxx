#include "itkTensorsToVTKCommand.h"

#include "itkTensorImageIO.h"
#include <itksys/SystemTools.hxx>

#include <vtkUnstructuredGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>

#include "GetPot.h"

namespace itk
{

  TensorsToVTKCommand::TensorsToVTKCommand()
  {
    m_ShortDescription = "Convert a tensor image (or a list of tensor images) into a vtkUnstructuredGrid structure";
    m_LongDescription = m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    m_LongDescription +="-i  [input tensor  image (or list of images (.list))]\n";    
    m_LongDescription +="-m  [optional mask image]\n";
    m_LongDescription +="-v  [output vectors instead of tensors (default 0)]\n";
    m_LongDescription +="-o  [output tensor mesh structure (vtk format)]\n";    
  }

  TensorsToVTKCommand::~TensorsToVTKCommand()
  {}

  int TensorsToVTKCommand::Execute (int narg, const char* arg[])
  {

    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char* input = cl.follow("nofile",2,"-I","-i");
    const char* output = cl.follow("nofile",2,"-O","-o");
    const char* mask = cl.follow("nofile",2,"-O","-o");
    const bool vectors = cl.follow(0,2,"-V","-v");
    
    std::cout << "Processing itk-vtk convert with following arguments: " << std::endl;
    std::cout << "input: " << input << std::endl;
    std::cout << "output: " << output << std::endl;
    if (cl.search(2,"-m","-M")) std::cout << "mask: " << mask << std::endl;
    std::cout << std::flush;

    // typedefs
    typedef double                                                         ScalarType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>                           TensorIOType;
    typedef TensorIOType::TensorImageType                                  TensorImageType;
    typedef TensorImageType::PixelType                                     TensorType;
    typedef itk::ImageRegionIterator<TensorImageType>                      TensorIteratorType;
    typedef itk::Matrix<ScalarType, 3, 3>                                  MatrixType;
    typedef itk::Vector<double, 3>                                         VectorType;
    typedef VectorType                                                     DisplacementType;
    typedef TensorImageType::PointType                                     PointType;
    typedef std::vector<TensorType> TensorArrayType;
    typedef std::vector<VectorType> VectorArrayType;
    typedef std::vector<PointType> PointArrayType;
    typedef itk::Image<unsigned short, 3> ImageType;
    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    
    ImageReaderType::Pointer maskreader = ImageReaderType::New();
    maskreader->SetFileName( mask );
    if (cl.search(2,"-m","-M"))
    {
      
      std::cout << "Reading: " << mask << std::flush;
      try
      {
	maskreader->Update();
      }
      catch( itk::ExceptionObject &e)
      {
	std::cerr << e;
	return -1;
      }
      std::cout << " Done." << std::endl;
    }
    
    ImageType::Pointer maskImage = maskreader->GetOutput();    
    
    std::vector<std::string> filelist;
    
    std::cout << "Reading input tensors: " << input << std::endl;
    
    std::cout<<"reading list : "<<input<<std::endl;
    std::string inputstring = input;
    std::string extension = itksys::SystemTools::GetFilenameLastExtension(inputstring).c_str();
    
    if (strcmp (extension.c_str(), ".mha") != 0)
    {
      std::ifstream inputliststream (input);
      if(inputliststream.fail())
      {
	std::cerr << "Unable to open file: " << input << std::endl;
	std::exit (EXIT_FAILURE);
      }
      unsigned int NumberOfImageFiles = 0;
      inputliststream >> NumberOfImageFiles;
      std::cout<<"encountered : "<<NumberOfImageFiles<<std::endl;
      
      std::string sline = "";
      itksys::SystemTools::GetLineFromStream(inputliststream, sline);
      
      for (unsigned int N=0; N<NumberOfImageFiles; N++)
      {
	std::string line = "";
	itksys::SystemTools::GetLineFromStream(inputliststream, line);
	filelist.push_back (line);
      }
    }
    else
      filelist.push_back (inputstring);
    
    TensorArrayType vecT;
    VectorArrayType vecV;
    PointArrayType vecP;
    
    for (unsigned int N=0; N<filelist.size(); N++)
    {
      std::string line = filelist[N];
      std::cout << "Reading tensors: " << line.c_str() << std::flush;
      TensorIOType::Pointer reader = TensorIOType::New();
      reader->SetFileName(line.c_str());
      bool success = false;
      
      try
      {
	reader->Read();
	success = true;
      }
      catch(itk::ExceptionObject &e)
      {
	std::cerr << "cannot read with TensorIO"<< std::endl;
      }
      
      if (!success)
      {
	vtkDataSetReader* vtkreader = vtkDataSetReader::New();
	vtkreader->SetFileName (line.c_str());
	vtkreader->Update();
	vtkPointSet* dataset = vtkPointSet::SafeDownCast (vtkreader->GetOutput());
	
	if (!dataset || !dataset->GetNumberOfPoints() || !dataset->GetPointData()->GetTensors())
	{
	  std::cerr << "cannot read file "<<line.c_str()<< std::endl;
	  std::cerr << "skipping line..."<< std::endl;
	  vtkreader->Delete();
	  
	  continue;
	}
	
	for (unsigned int i=0; i<dataset->GetNumberOfPoints(); i++)
	{
	  PointType x;
	  for (unsigned int j=0; j<3; j++)
	    x[j] = dataset->GetPoint (i)[j];
	  TensorType T;
	  VectorType V;
	  double* vals = dataset->GetPointData()->GetTensors()->GetTuple (i);
	  
	  T[0] = vals[0];
	  T[1] = vals[1];
	  T[2] = vals[4];
	  T[3] = vals[2];
	  T[4] = vals[5];
	  T[5] = vals[8];
	  
	  V = T.GetEigenvector (2);
	  V.Normalize();
	  
	  vecT.push_back(T);
	  vecV.push_back (V);
	  
	  vecP.push_back(x);
	}
	
	vtkreader->Delete();
	
      }
      
      TensorImageType::Pointer tensorimage = reader->GetOutput();
      TensorIteratorType it(tensorimage, tensorimage->GetLargestPossibleRegion());
      
      itk::ImageRegionConstIterator<ImageType> itMask (maskImage, maskImage->GetLargestPossibleRegion() );
      
      while( !it.IsAtEnd() )
      {
	PointType x;
	tensorimage->TransformIndexToPhysicalPoint(it.GetIndex(), x);
	
	TensorType T = it.Get();
	VectorType V;
	
	bool isinside = 1;
	
	if (cl.search(2,"-m","-M"))
	{
	  ImageType::IndexType index;
	  isinside = maskImage->TransformPhysicalPointToIndex (x, index);
	  if (isinside)
	  {
	    itMask.SetIndex (index);
	    isinside = itMask.Value() > 0;
	  }
	}
	
	if (isinside && (T.GetTrace() > 0.0) && !T.HasNans())
	{	
	  T = T.ApplyMatrix(tensorimage->GetDirection());
	  V = T.GetEigenvector (2);
	  V.Normalize();
	  
	  vecT.push_back(T);
	  vecV.push_back(V);
	  vecP.push_back(x);
	}
	++it;
      }
      std::cout<<" done."<<std::endl; 
    }
    
    vtkUnstructuredGrid* crossvalidation = vtkUnstructuredGrid::New();
    vtkDoubleArray* data1 = vtkDoubleArray::New();
    vtkDoubleArray* data2 = vtkDoubleArray::New();
    
    vtkPoints* points = vtkPoints::New();
    points->SetNumberOfPoints (vecP.size());
    
    data1->SetNumberOfComponents (9);
    data1->SetNumberOfTuples (vecP.size());
    data2->SetNumberOfComponents (3);
    data2->SetNumberOfTuples (vecP.size());
    
    for (unsigned int i=0; i<vecP.size(); i++)
    {
      double pt[3] = {vecP[i][0],
		      vecP[i][1],
		      vecP[i][2]};
      
      points->SetPoint (i, pt);
      
      double vals1[9];
      vals1[0] = vecT[i][0];
      vals1[1] = vecT[i][1];
      vals1[2] = vecT[i][3];
      vals1[3] = vecT[i][1];
      vals1[4] = vecT[i][2];
      vals1[5] = vecT[i][4];
      vals1[6] = vecT[i][3];
      vals1[7] = vecT[i][4];
      vals1[8] = vecT[i][5];
      
      double vals2[3];      
      vals2[0] = vecV[i][0];
      vals2[1] = vecV[i][1];
      vals2[2] = vecV[i][2];
      
      data1->SetTuple (i, vals1);
      data2->SetTuple (i, vals2);
    }
    
    crossvalidation->SetPoints (points);
    if (vectors)
      crossvalidation->GetPointData()->SetVectors (data2);
    else
      crossvalidation->GetPointData()->SetTensors (data1);
    
    vtkDataSetWriter* writer = vtkDataSetWriter::New();
    writer->SetFileName (output);
    writer->SetInput (crossvalidation);
    writer->Update();
    
    writer->Delete();
    points->Delete();
    data1->Delete();
    data2->Delete();
    
    crossvalidation->Delete();
    
    std::cout<<" done."<<std::endl;    
    
    return EXIT_SUCCESS;
  }
  
}
