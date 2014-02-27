/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkApplyTransformToImage.cxx 1 2010-05-21 14:00:33Z nt08 $
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

/**
 * 
 * \file itkApplyTransformToImage.cxx
 *
 * \brief  Apply a simple transform to an image, in ITK format. 
 * 
 * This tool takes an image and a transform in ITK format as inputs.
 * The transform is applied to the "header" of the image, meaning
 * it will only transform the origin and 3x3 direction of the input
 * image. Therefore, the input transform has to be a derivative of
 * a MatrixOffsetTransformBase, i.e. an affine transform
 *
 * \note This tool does *NOT* check if the input transform matrix
 * is rigid. If it is not the case the output image could have direction
 * vectors that are *NOT* orthogonal to eachother. 
 * 
 * \author Nicolas Toussaint, King's College London, nicolas.toussaint@kcl.ac.uk
 * 
 * \url http://www.kcl.ac.uk/schools/medicine/research/imaging/
 *  
 */
#include <itkApplyTransformToMeshCommand.h>


#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkPointSet.h>

#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include "itkTranslationTransform.h"
#include <itkMatrixOffsetTransformBase.h>

#include "GetPot.h"


namespace itk
{

  ApplyTransformToMeshCommand::ApplyTransformToMeshCommand()
  {
    m_ShortDescription = "Apply an affine transformation to a vtk mesh";
    m_LongDescription = m_ShortDescription;
    m_LongDescription = "\n\nUsage:\n";
    m_LongDescription +="-i [input image]\n";
    m_LongDescription +="-t [transform]\n";
    m_LongDescription +="-o [Output file name]\n";
  }

  ApplyTransformToMeshCommand::~ApplyTransformToMeshCommand()
  {}

  int ApplyTransformToMeshCommand::Execute (int narg, const char* arg[])
  {


    GetPot cl(narg, const_cast<char**>(arg)); // argument parser
    if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    {
      std::cout << std::endl << this->GetLongDescription() << std::endl;
      return -1;
    }
    
    const bool IsInputPresent    = cl.search(2,"-i","-I");
    const bool IsOutputPresent   = cl.search(2,"-o","-O");
  
    if(!IsInputPresent || !IsOutputPresent )
    {
      std::cerr << "Error: Input and (or) output not set." << std::endl;
      exit (EXIT_FAILURE);
    }


    const char* fileOut  = cl.follow("NoFile",2,"-o","-O");
    const char* fileIn  = cl.follow("NoFile",2,"-i","-I");
    const char* fileTr  = cl.follow("NoFile",2,"-t","-T");

    
    typedef double                                           ScalarType;  
    typedef itk::MatrixOffsetTransformBase<ScalarType, 3, 3> LinearTransformType;
    typedef itk::TranslationTransform<ScalarType, 3>         TranslationTransformType;
    typedef itk::Transform<double, 3, 3>                     TransformType;
    typedef itk::TransformFileReader                         TransformReaderType;

    itk::TransformFactory< LinearTransformType >::RegisterTransform ();
    itk::TransformFactory< TranslationTransformType >::RegisterTransform ();

    TransformReaderType::Pointer transforeader = TransformReaderType::New();
    transforeader->SetFileName ( fileTr );
    try
    {
      transforeader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return EXIT_FAILURE;
    }
  
    LinearTransformType::Pointer transform1 = dynamic_cast<LinearTransformType*>( transforeader->GetTransformList()->front().GetPointer() );
    
    if (!transform1)
    {
      std::cerr << "The transformation written in "<<fileTr<<" does not derive from a MatrixOffsetTransformBase, "
		<<"which is mandatory in this executable to be able to only change the header of the file"<<std::endl;
      return EXIT_FAILURE;
    }
    
    vtkDataSetReader* reader = vtkDataSetReader::New();
    reader->SetFileName (fileIn);
    reader->Update();
    vtkPointSet* pointset   = vtkPointSet::SafeDownCast (reader->GetOutput());
    vtkPoints* inputpoints  = pointset->GetPoints();
    vtkPoints* outputpoints = vtkPoints::New();
    vtkMatrix4x4* matrix = vtkMatrix4x4::New(); // identity by default
    matrix->Identity();
    
    // somehow we do not have to correct for the orientation
    // in the specific case demanded
    
    for (unsigned int x=0; x<4; x++)
      for (unsigned int y=0; y<4; y++)
	matrix->SetElement (x,y,transform1->GetMatrix()[x][y]);
    std::cout<<"the matrix is as follows :"<<std::endl;
    std::cout<<(*matrix)<<std::endl;
    
    vtkTransform* transform = vtkTransform::New();
    transform->SetMatrix (matrix);
    transform->TransformPoints (inputpoints, outputpoints);
    pointset->SetPoints (outputpoints);
    vtkDataSetWriter* writer = vtkDataSetWriter::New();
    // writer->SetFileTypeToBinary();
    writer->SetFileName (fileOut);
    writer->SetInputData (pointset);
    writer->Update();
    reader->Delete();
    outputpoints->Delete();
    transform->Delete();
    matrix->Delete();
    writer->Delete();

    std::cout << " Done." << std::endl;  
  
    return EXIT_SUCCESS;

  }
}
