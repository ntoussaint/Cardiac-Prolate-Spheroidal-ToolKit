/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkExtrapolateTensorField.cxx 1 2010-05-21 14:00:33Z nt08 $
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
#include "itkCumulateProlateSpheroidalDataSetsCommand.h"
#include "itkTensorMeshIO.h"
#include <itksys/SystemTools.hxx>

#include "GetPot.h"


typedef double                                                         ScalarType;
typedef itk::TensorMeshIO <ScalarType, 3, 3>                           TensorMeshIOType;
typedef TensorMeshIOType::TensorMeshType                               MeshType;


template<typename MeshType>
void EvaluateXi1Range(typename MeshType::Pointer mesh, double* range)
{
  double min = range[0];
  double max = range[1];
  typename MeshType::PointType p;
  for (unsigned int i=0; i<mesh->GetNumberOfPoints(); i++)
  {
    mesh->GetPoint (i, &p);
    min = std::min (p[0], min);
    max = std::max (p[0], max);    
  }
  range[0] = min;
  range[1] = max;
  
}


template<typename MeshType>
void NormalizeXi1(typename MeshType::Pointer mesh, double* outrange)
{
  double inrange[2] = {+100, -100};
  EvaluateXi1Range<MeshType> (mesh, inrange);
  
  double inmin = inrange[0];
  double inmax = inrange[1];
  double outmin = outrange[0];
  double outmax = outrange[1];
  typename MeshType::PointType p;
  
  std::cout<<"1st component normalization:"<<std::endl
	   <<"from : "<<inmin<<" and "<<inmax<<std::endl
	   <<"  to : "<<outmin<<" and "<<outmax<<std::endl;
  
  for (unsigned int i=0; i<mesh->GetNumberOfPoints(); i++)
  {
    mesh->GetPoint (i, &p);
    p[0] = outmin + ( (outmax - outmin) / (inmax - inmin) ) * (p[0] - inmin);
    mesh->SetPoint (i, p);
  }
}


template<typename MeshType>
void AppendMesh(typename MeshType::Pointer in, typename MeshType::Pointer out)
{
  unsigned int outindex = out->GetNumberOfPoints();
  typename MeshType::PointType p;
  typename MeshType::PixelType t;
  for (unsigned int inindex = 0; inindex<in->GetNumberOfPoints(); inindex++)
  {
    in->GetPoint (inindex, &p);
    in->GetPointData (inindex, &t);
    out->GetPoints()->InsertElement (outindex, p);
    out->GetPointData()->InsertElement (outindex, t);
    outindex++;
  }
}


namespace itk
{

  CumulateProlateSpheroidalDataSetsCommand::CumulateProlateSpheroidalDataSetsCommand()
  {
    m_ShortDescription = "Cumulate several datasets in one, and normalizing 1st component";
    m_LongDescription += m_ShortDescription;
    m_LongDescription += "\nThe first component is normalized against the first mesh in the list";
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-i    [input list of tensor meshes]\n";
    m_LongDescription += "-o    [output tensor cumulated mesh (use 1st mesh for back transformation) \n";
  }
  
  CumulateProlateSpheroidalDataSetsCommand::~CumulateProlateSpheroidalDataSetsCommand()
  {}

  int CumulateProlateSpheroidalDataSetsCommand::Execute (int narg, const char* arg[])
  {

    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
  
    const char*  inputfile                    = cl.follow("NoFile",2,"-i","-I");
    const char*  outputfile                    = cl.follow("NoFile",2,"-o","-O");
    
    std::cout << "Processing atlas creation: " << std::endl;
    std::cout << std::flush;
  
    // read the input tensors and put tham into a vtkUnstructuredGrid
  
    std::cout << "Reading input tensors list: " << inputfile << std::endl;
    std::vector<std::string> filelist;
    std::string inputstring = inputfile;
    std::ifstream inputliststream (inputfile);
    if(inputliststream.fail())
    {
      std::cerr << "Unable to open file: " << inputfile << std::endl;
      std::exit (EXIT_FAILURE);
    }
    unsigned int NumberOfFiles = 0;
    inputliststream >> NumberOfFiles;
    std::cout<<"number of files : "<<NumberOfFiles<<std::endl;
    std::string sline = "";
    itksys::SystemTools::GetLineFromStream(inputliststream, sline);
    
    for (unsigned int N=0; N<NumberOfFiles; N++)
    {
      std::string line = "";
      itksys::SystemTools::GetLineFromStream(inputliststream, line);
      filelist.push_back (line);
    }
    
    MeshType::Pointer                        data = MeshType::New();
    MeshType::PointsContainer::Pointer     points = MeshType::PointsContainer::New();
    MeshType::PointDataContainer::Pointer tensors = MeshType::PointDataContainer::New();
    data->SetPoints (points);
    data->SetPointData (tensors);
    
    double xi1range[2] = {+100, -100};
    
    for (unsigned int i=0; i<filelist.size(); i++)
    {
      std::cout<<"appending : "<<filelist[i].c_str()<<" ..."<<std::flush;
      TensorMeshIOType::Pointer reader = TensorMeshIOType::New();
      reader->SetFileName (filelist[i].c_str());
      reader->Read();
      MeshType::Pointer mesh = reader->GetOutput();
      if (i == 0)
	EvaluateXi1Range<MeshType> (mesh, xi1range);
      else
	NormalizeXi1<MeshType> (mesh, xi1range);
      
      AppendMesh<MeshType> (mesh, data);
      std::cout<<" Done."<<std::endl;
    }

    TensorMeshIOType::Pointer writer = TensorMeshIOType::New();
    writer->SetInput (data);
    writer->SetFileName (outputfile);
    writer->Write();
      
    return EXIT_SUCCESS;
  }

}
