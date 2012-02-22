/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkCreateDomain.cxx 1 2010-05-21 14:00:33Z nt08 $
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
#include <itkFindClosestProlateTransformCommand.h>

/**

 * 
 *
 * \author Nicolas Toussaint, King's College London
 */
#include <vtkMath.h>
#include <vtkPCAAnalysisFilter.h>
#include <vtkMatrix3x3.h>
#include <vtkMatrix4x4.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkMath.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkTransformFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include <itkMatrix.h>

#include "itkProlateSpheroidalTransform.h"
#include "itkTensor.h"

#include <sstream>
#include <fstream>
#include <cstdlib>
#include "GetPot.h"

namespace itk
{

  FindClosestProlateTransformCommand::FindClosestProlateTransformCommand()
  {
    m_ShortDescription = "Find a prolate spheroid approximating an LV segmentation mesh";
    m_LongDescription += m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-i    [input midwall (or total LV) mesh file]\n";
    m_LongDescription += "-p    [percentage of basal ventricle above 90 deg. (default: 0.8)]\n";
    m_LongDescription += "-o    [output prolate transform (default prolatetransform.tr)]\n";
    m_LongDescription += "-v    [verbose]\n";
  }

  FindClosestProlateTransformCommand::~FindClosestProlateTransformCommand()
  {}
  
  int FindClosestProlateTransformCommand::Execute (int narg, const char* arg[])
  {
    
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const char*  inputfile  = cl.follow("nofile",2,"-i","-I");
    const char*  outputfile = cl.follow("prolatetransform.tr",2,"-o","-O");
    const double percentile = cl.follow(0.8,2,"-p","-P");
    const bool   verbose    = cl.follow(false,2,"-v","-V");
  
    std::cout << "Processing prolate search with following arguments: " << std::endl;
    std::cout << "inputfile: " << inputfile << std::endl;
    std::cout << "percentile: " << percentile << std::endl;
    std::cout << "output: " << outputfile << std::endl;
    std::cout << "verbose: " << verbose << std::endl;
    std::cout << std::flush;
  
    typedef double ScalarType;
    typedef itk::Image<ScalarType, 3> ImageType;
    typedef itk::ImageFileReader <ImageType> ImageReaderType; 
    typedef itk::ImageFileWriter <ImageType> ImageWriterType;
    typedef itk::ImageFileWriter <ImageType> DoubleImageWriterType;
    typedef itk::ProlateSpheroidalTransform<ScalarType> TransformType;
    typedef TransformType::InputPointType PointType;
  
    vtkDataSetReader* reader = vtkDataSetReader::New();
    reader->SetFileName (inputfile);
    reader->Update();

    vtkPointSet* midwall = vtkPointSet::SafeDownCast (reader->GetOutput());
    if (!midwall)
    {
      std::cerr<<"input not a vtk pointset class, exiting"<<std::endl;
      reader->Delete();
      std::exit (EXIT_FAILURE);
    }
    int numberofpoints = midwall->GetNumberOfPoints();

    /// compute the ventricle's center of mass
    itk::Point<double,3> point;
    itk::Point<double,3> mean; mean[0] = mean[1] = mean[2] = 0.0;
    for (int i=0; i<numberofpoints; i++)
    {
      midwall->GetPoint (i, point.GetDataPointer());
      mean += point.GetVectorFromOrigin();
    }
    for (int i=0; i<3; i++) mean[i] /= (double)numberofpoints;

    /// compute the covariance matrix of all ventricle's points
    vnl_matrix_fixed<double,3,3> matrix (0.0);  
    for (int i=0; i<numberofpoints; i++)
    {
      midwall->GetPoint (i, point.GetDataPointer());
      for (unsigned int j=0; j<3; j++)
	for (unsigned int k=0; k<3; k++)
	  matrix[j][k] += (point[j] - mean[j])*(point[k] - mean[k]);
    }
    matrix /= (double)numberofpoints;

    itk::Tensor<double,3> covariance;
    covariance.SetVnlMatrix (matrix);

    /// compute the long axis, first eigen vector of the covariance matrix
    itk::Vector<double,3> longaxis = covariance.GetEigenvector (2);
    longaxis.Normalize();

    /// compute coordinates of apex and base
    double projectionrange[2] = {VTK_DOUBLE_MAX, VTK_DOUBLE_MIN};
    itk::Point<double,3> apex, base;
    apex[0] = apex[1] = apex[2] = 0.0;
    base[0] = base[1] = base[2] = 0.0;
  
    itk::Vector<double,3> vec;
    for (int i=0; i<numberofpoints; i++)
    {
      midwall->GetPoint (i, point.GetDataPointer());
      vec = (point - mean);
      double projection = vec * longaxis;
      if (projection > projectionrange[1])
      {
	base = mean + ( percentile * projection ) * longaxis;
	projectionrange[1] = projection;
      }
      if (projection < projectionrange[0])
      {
	apex = mean + ( percentile * projection ) * longaxis;
	projectionrange[0] = projection;
      }
    }

    /// compute coordinates of the ventricle's radius at base --> septum
    double middle = (projectionrange[0] + projectionrange[1]) / 2.0;
    double radius = 0;
    int midnumberofpoints = 0;
    for (int i=0; i<numberofpoints; i++)
    {
      midwall->GetPoint (i, point.GetDataPointer());
      vec = point - mean;
      double projection = vec * longaxis;
      if (projection < middle)
	continue;
      radius += std::sqrt ( vec.GetSquaredNorm() - projection * projection );
      midnumberofpoints ++;
    }
    radius /= (double)midnumberofpoints;
    itk::Point<double,3> septum = base + radius * covariance.GetEigenvector (1);

    /// Fill the Prolate Spheroidal Transform with information
    TransformType::ParametersType parameters;
    parameters.SetSize (9);
    parameters[0] = base[0];   parameters[1] = base[1];   parameters[2] = base[2];
    parameters[3] = apex[0];   parameters[4] = apex[1];   parameters[5] = apex[2]; 
    parameters[6] = septum[0]; parameters[7] = septum[1]; parameters[8] = septum[2];
    TransformType::Pointer transform = TransformType::New();
    transform->SetParameters (parameters);
    
    /// Write the transform
    itk::TransformFactory<TransformType>::RegisterTransform ();
    itk::TransformFileWriter::Pointer transformwriter = itk::TransformFileWriter::New();
    transformwriter->SetFileName(outputfile);
    transformwriter->SetInput (transform);
    transformwriter->Update();

    if (verbose)
    {
      std::cout<<"base, apex, septum"<<std::endl;
      std::cout<<base<<std::endl;
      std::cout<<apex<<std::endl;
      std::cout<<septum<<std::endl;
    }

    reader->Delete();
  
    return EXIT_SUCCESS;
  }

}
