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
#include "itkCreateProlateDomainCommand.h"

/**

 * 
 *
 * \author Nicolas Toussaint, King's College London
 */
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <vtkPointLocator.h>
#include <vtkMath.h>

#include <vtkMatrixToLinearTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkUnstructuredGrid.h>

#include "vtkDataSetReader.h"
#include "vtkDataSetWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "vtkPolyData.h"
#include "vtkImageData.h"
#include <vtkImageThreshold.h>
#include <sstream>
#include <fstream>
#include <cstdlib>

#include "itkProlateSpheroidalTransform.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include <itkMultiplyImageFilter.h>

#include "GetPot.h"

typedef unsigned short BinaryScalarType;
typedef itk::Image<BinaryScalarType, 3> BinaryImageType;
typedef double ScalarType;
typedef itk::Image<ScalarType, 3> ImageType;
typedef itk::ImageFileReader <BinaryImageType> ImageReaderType; 
typedef itk::ImageFileWriter <BinaryImageType> ImageWriterType;
typedef itk::ImageFileWriter <ImageType> DoubleImageWriterType;
typedef itk::ProlateSpheroidalTransform<ScalarType> TransformType;
typedef TransformType::InputPointType PointType;
typedef itk::BinaryThresholdImageFilter <ImageType, BinaryImageType> ThresholdType;
typedef itk::MultiplyImageFilter <BinaryImageType, BinaryImageType> MultiplyType;

void CreateGrid (TransformType::Pointer transform,
		 double mu1, double mu2,
		 double nu1, double nu2,
		 unsigned int throughwalldivisions,
		 unsigned int longdivisions,
		 unsigned int circumdivisions,
		 vtkUnstructuredGrid* grid,
		 bool stay_in_prolate = 0,
		 bool imitatecircular = 0)
{

  unsigned int N[3], N1[3], N2[3];
  N[0] = throughwalldivisions + 1;
  N[1] = longdivisions + 1;
  N[2] = circumdivisions + 1;
  
  // unsigned int start_long_axis = itk::Math::Round ((double)N[1] / 6.0);
  // unsigned int stop_long_axis = itk::Math::Round ((double)N[1] / 2.5);
  // unsigned int start_circum = 0;
  // unsigned int stop_circum = itk::Math::Round ((double)N[2] / 3.0);
  
  // unsigned int start_long_axis = 0;
  // unsigned int stop_long_axis = N[1];
  // unsigned int start_circum = 0;
  // unsigned int stop_circum = N[2];

  unsigned int start_long_axis = 0;
  unsigned int stop_long_axis = N[1];
  // unsigned int start_circum = itk::Math::Round ((double)N[2] * 0.0 / 50.0);
  // unsigned int stop_circum = itk::Math::Round ((double)N[2] * 12.5 / 50.0);
  unsigned int start_circum = 0;
  unsigned int stop_circum = N[2];
    
  N1[0] = 0;
  N1[1] = start_long_axis;
  N1[2] = start_circum;
  
  N2[0] = N[0];
  N2[1] = stop_long_axis;
  N2[2] = stop_circum;
  
  double bounds[3][2];
  bounds[0][0] = mu1; bounds[0][1] = mu2;
  bounds[1][0] = nu1; bounds[1][1] = nu2;
  bounds[2][0] = 0.0; bounds[2][1] = 2 * vnl_math::pi;
  
  vtkPoints* outputpoints = vtkPoints::New();
  
  grid->SetPoints (outputpoints);
  grid->Allocate();

  double zeta1factor = 40;
  double zeta2factor = 4;

  for (unsigned int i=N1[0]; i<N2[0]; i++)
    for (unsigned int j=N1[1]; j<N2[1]; j++)
      for (unsigned int k=N1[2]; k<N2[2]; k++)
      {
	vtkIdType cell[2];
	
	PointType point1, point2, pt;
	point1[0] = bounds[0][0] + (double)i * ( (bounds[0][1] - bounds[0][0]) / (double)(N[0] - 1) );
	point1[1] = bounds[1][0] + (double)j * ( (bounds[1][1] - bounds[1][0]) / (double)(N[1] - 1) );
	point1[2] = bounds[2][0] + (double)k * ( (bounds[2][1] - bounds[2][0]) / (double)(N[2] - 1) );
	if (point1[2] > 2.0 * vnl_math::pi) point1[2] =  2.0 * vnl_math::pi - vcl_numeric_limits<double>::epsilon();
	
	// pt = transform->TransformPoint (point1);
	if (!stay_in_prolate)
	  pt = transform->TransformPoint (point1);
	else
	{
	  pt[0] = zeta1factor*point1[0];
	  if (imitatecircular)
	  {
	    pt[1] = zeta2factor*point1[1]*std::cos (point1[2]);
	    pt[2] = zeta2factor*point1[1]*std::sin (point1[2]);
	  }
	  else
	  {
	    pt[1] = zeta2factor*point1[1];
	    pt[2] = point1[2];	    
	  }
	}
	
	cell[0] = outputpoints->InsertNextPoint (pt.GetDataPointer());
	
	if (i <N2[0]-1)
	{
	  point2[0] = bounds[0][0] + (double)(i+1) * ( (bounds[0][1] - bounds[0][0]) / (double)(N[0] - 1) );
	  point2[1] = bounds[1][0] + (double)j * ( (bounds[1][1] - bounds[1][0]) / (double)(N[1] - 1) );
	  point2[2] = bounds[2][0] + (double)k * ( (bounds[2][1] - bounds[2][0]) / (double)(N[2] - 1) );
	  if (point2[2] > 2.0 * vnl_math::pi) point2[2] =  2.0 * vnl_math::pi - vcl_numeric_limits<double>::epsilon();

	  // pt = transform->TransformPoint (point2);
	  if (!stay_in_prolate)
	    pt = transform->TransformPoint (point2);
	  else
	  {
	    pt[0] = zeta1factor*point2[0];
	    if (imitatecircular)
	    {
	      pt[1] = zeta2factor*point2[1]*std::cos (point2[2]);
	      pt[2] = zeta2factor*point2[1]*std::sin (point2[2]);
	    }
	    else
	    {
	      pt[1] = zeta2factor*point2[1];
	      pt[2] = point2[2];    
	    }	      
	  }
	  
	  cell[1] = outputpoints->InsertNextPoint (pt.GetDataPointer());
	  grid->InsertNextCell (VTK_LINE, 2, cell);
	}
	if (j <N2[1]-1)
	{
	  point2[0] = bounds[0][0] + (double)i * ( (bounds[0][1] - bounds[0][0]) / (double)(N[0] - 1) );
	  point2[1] = bounds[1][0] + (double)(j+1) * ( (bounds[1][1] - bounds[1][0]) / (double)(N[1] - 1) );
	  point2[2] = bounds[2][0] + (double)k * ( (bounds[2][1] - bounds[2][0]) / (double)(N[2] - 1) );
	  if (point2[2] > 2.0 * vnl_math::pi) point2[2] =  2.0 * vnl_math::pi - vcl_numeric_limits<double>::epsilon();

	  // pt = transform->TransformPoint (point2);
	  if (!stay_in_prolate)
	    pt = transform->TransformPoint (point2);
	  else
	  {
	    pt[0] = zeta1factor*point2[0];
	    if (imitatecircular)
	    {
	      pt[1] = zeta2factor*point2[1]*std::cos (point2[2]);
	      pt[2] = zeta2factor*point2[1]*std::sin (point2[2]);
	    }
	    else
	    {
	      pt[1] = zeta2factor*point2[1];
	      pt[2] = point2[2];
	    }
	  }
	  
	  cell[1] = outputpoints->InsertNextPoint (pt.GetDataPointer());
	  grid->InsertNextCell (VTK_LINE, 2, cell);
	}
	if (k <N2[2]-1)
	{
	  point2[0] = bounds[0][0] + (double)i * ( (bounds[0][1] - bounds[0][0]) / (double)(N[0] - 1) );
	  point2[1] = bounds[1][0] + (double)j * ( (bounds[1][1] - bounds[1][0]) / (double)(N[1] - 1) );
	  point2[2] = bounds[2][0] + (double)(k+1) * ( (bounds[2][1] - bounds[2][0]) / (double)(N[2] - 1) );
	  if (point2[2] > 2.0 * vnl_math::pi) point2[2] =  2.0 * vnl_math::pi - vcl_numeric_limits<double>::epsilon();

	  // pt = transform->TransformPoint (point2);
	  if (!stay_in_prolate)
	    pt = transform->TransformPoint (point2);
	  else
	  {
	    pt[0] = zeta1factor*point2[0];
	    if (imitatecircular)
	    {
	      pt[1] = zeta2factor*point2[1]*std::cos (point2[2]);
	      pt[2] = zeta2factor*point2[1]*std::sin (point2[2]);
	    }
	    else
	    {
	      pt[1] = zeta2factor*point2[1];
	      pt[2] = point2[2];
	    }
	  }
	  
	  cell[1] = outputpoints->InsertNextPoint (pt.GetDataPointer());
	  grid->InsertNextCell (VTK_LINE, 2, cell);
	}
      }

  outputpoints->Delete();
  
}



namespace itk
{

  CreateProlateDomainCommand::CreateProlateDomainCommand()
  {
    m_ShortDescription = "Create an ellipsoidal domain image out of a prolate sph. transform";
    m_LongDescription += m_ShortDescription;
    m_LongDescription += "\n\nUsage:\n";
    
    m_LongDescription += "-i    [input  geometry (i.e. anatomy/domain)]\n";
    m_LongDescription += "-t    [prolate transform]\n";
    m_LongDescription += "-p    [phase (0: diastole / 1: systole)]\n";
    m_LongDescription += "-o    [output ellipsoidal image]\n";
    m_LongDescription += "-w    [wall thickness in mm (default: diastole 12 / systole 17)]\n";
    m_LongDescription += "-a    [basal angle in deg. (default: 95 deg.)]\n";
  }

  CreateProlateDomainCommand::~CreateProlateDomainCommand()
  {}
  
  
  int CreateProlateDomainCommand::Execute (int narg, const char* arg[])
  {
    
    GetPot   cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") ) 
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    std::string domainfile = cl.follow("nofile",2,"-i","-I");
    std::string transformfile = cl.follow("nofile",2,"-t","-T");
    std::string ellipsoiddomainfile = cl.follow("nofile",2,"-o","-O");
    const bool phase = cl.follow(0,2,"-p","-P");
    
    double thickness = 0.0;
    switch (phase)
    {
	case 1: // systole
	  thickness = 17.0;
	  break;
	case 0: // diastole
	default:
	  thickness = 12.0;
    }

    if (cl.search(2,"-w","-W"))
      thickness = cl.follow(12.0,2,"-w","-W");
    
    double maxangle = 95.0;
    if (cl.search(2,"-a","-A"))
      maxangle = cl.follow(95.0,2,"-a","-A");    

    std::cout<<"found thickness of "<<thickness<<std::endl;
    std::cout<<"found maxangle of "<<maxangle<<std::endl;
    
    
    // read the input image
    ImageReaderType::Pointer imReader = ImageReaderType::New();
    imReader->SetFileName(domainfile.c_str());
    std::cout << "Reading: " << domainfile.c_str() << std::flush;
    try
    {
      imReader->Update();
    }
    catch(itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << " Done." << std::endl;
    BinaryImageType::Pointer domain = imReader->GetOutput();
  
    std::cout<<"reading transform file "<<transformfile.c_str()<<std::endl;
    TransformType::Pointer transform = TransformType::New();
    itk::TransformFactory<TransformType>::RegisterTransform ();
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName(transformfile.c_str());
    transformreader->Update();
    transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    std::cout<<"done."<<std::endl;

    double mu1 = TransformType::asinh ((transform->GetLambda2() - thickness/2.0) / transform->GetSemiFociDistance());
    double mu2 = TransformType::asinh ((transform->GetLambda2() + thickness/2.0) / transform->GetSemiFociDistance());
    double nu1 = 0.0;
    double nu2 = maxangle * vnl_math::pi / 180.0;
  
    std::cout<<"mu1 and mu2 "<<mu1<<" and "<<mu2<<std::endl;
    std::cout<<"nu1 and nu2 "<<nu1<<" and "<<nu2<<std::endl;

    // unsigned int circumdivisions = 14;
    unsigned int circumdivisions = 50;
    unsigned int longdivisions = 6;
    unsigned int throughwalldivisions = 3;
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    TransformType::Pointer inversetransform = TransformType::New();
    transform->GetInverse (inversetransform);
  
    CreateGrid (inversetransform,
    		mu1, mu2,
    		nu1, nu2,
    		throughwalldivisions,
    		longdivisions,
    		circumdivisions,
    		grid);

    vtkDataSetWriter* gridwriter = vtkDataSetWriter::New();
    gridwriter->SetInputData (grid);
    gridwriter->SetFileName ("grid.vtk");
    gridwriter->Update();
    grid->Delete();

    vtkUnstructuredGrid* grid2 = vtkUnstructuredGrid::New();
  
    CreateGrid (inversetransform,
    		mu1, mu2,
    		nu1, nu2,
    		throughwalldivisions,
    		longdivisions,
    		circumdivisions,
    		grid2, 1);

    gridwriter->SetInputData (grid2);
    gridwriter->SetFileName ("grid-prolate.vtk");
    gridwriter->Update();
  
    gridwriter->Delete();
    grid2->Delete();

    ImageType::Pointer zeta1 = ImageType::New();
    zeta1->SetRegions (domain->GetLargestPossibleRegion());
    zeta1->SetSpacing(domain->GetSpacing());
    zeta1->SetOrigin(domain->GetOrigin());
    zeta1->SetDirection(domain->GetDirection());
    zeta1->Allocate();
    ImageType::Pointer zeta2 = ImageType::New();
    zeta2->SetRegions (domain->GetLargestPossibleRegion());
    zeta2->SetSpacing(domain->GetSpacing());
    zeta2->SetOrigin(domain->GetOrigin());
    zeta2->SetDirection(domain->GetDirection());
    zeta2->Allocate();
  
    itk::ImageRegionIterator<ImageType> itzeta1(zeta1, zeta1->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType> itzeta2(zeta2, zeta2->GetLargestPossibleRegion());
    PointType in;
    PointType out;
    ImageType::PointType x;
  
    while( !itzeta1.IsAtEnd() )
    {
      zeta1->TransformIndexToPhysicalPoint(itzeta1.GetIndex(), x);      
      for (unsigned int i=0; i<3; i++)
	in[i] = x[i];
      out = transform->TransformPoint (in);
      itzeta1.Set (static_cast<ScalarType>(out[0]));
      itzeta2.Set (static_cast<ScalarType>(out[1]));
      ++itzeta1;
      ++itzeta2;
    }
  
    // DoubleImageWriterType::Pointer zetawriter = DoubleImageWriterType::New();
    // zetawriter->SetInput (zeta1);
    // zetawriter->SetFileName ("zeta1.mha");
    // zetawriter->Update();
  
  
    ThresholdType::Pointer threshold1 = ThresholdType::New();
    threshold1->SetInput(zeta1);
    threshold1->SetInsideValue (static_cast<BinaryScalarType>(1.0));
    threshold1->SetOutsideValue (static_cast<BinaryScalarType>(0.0));
    threshold1->SetUpperThreshold (mu2);
    threshold1->SetLowerThreshold (mu1);
    threshold1->Update();
    ThresholdType::Pointer threshold2 = ThresholdType::New();
    threshold2->SetInput(zeta2);
    threshold2->SetInsideValue (static_cast<BinaryScalarType>(1.0));
    threshold2->SetOutsideValue (static_cast<BinaryScalarType>(0.0));
    threshold2->SetUpperThreshold (nu2);
    threshold2->SetLowerThreshold (nu1);
    threshold2->Update();

    MultiplyType::Pointer multiplier = MultiplyType::New();
    multiplier->SetInput (0, threshold1->GetOutput());
    multiplier->SetInput (1, threshold2->GetOutput());
    multiplier->Update();
  
    ImageWriterType::Pointer ellipsoidwriter = ImageWriterType::New();
    ellipsoidwriter->SetInput (multiplier->GetOutput());
    ellipsoidwriter->SetFileName (ellipsoiddomainfile.c_str());
    ellipsoidwriter->Update();
    
    return EXIT_SUCCESS;
  
  
  }

}

