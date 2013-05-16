/*=========================================================================

  Program:   ImagingSciences
  Module:    $Id: itkApplyTransformToTensors.cxx 1 2010-05-21 14:00:33Z nt08 $
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
 * \file itkApplyTransformToTensors.cxx
 *
 * \brief  Apply a simple transform to a tensor image, in ITK format. 
 * 
 * 
 * This tool takes an tensor image and a transform in ITK format as inputs.
 * The transform is applied to the "header" of the tensor image, meaning
 * it will only transform the origin and 3x3 direction of the input
 * image. Therefore, the input transform has to be a derivative of
 * a MatrixOffsetTransformBase, i.e. an affine transform
 *
 * \note This tool does *NOT* check if the input transform matrix
 * is rigid. If it is not the case the output tensor image could have direction
 * vectors that are *NOT* orthogonal to eachother. 
 *
 *  
 * \author Nicolas Toussaint, King's College London, nicolas.toussaint@kcl.ac.uk
 * 
 * \url http://www.kcl.ac.uk/schools/medicine/research/imaging/
 *  
 */
#include <itkApplyTransformToTensorImageCommand.h>

#include <itkImage.h>
#include <cstdlib>
#include "itkTensorImageIO.h"
#include <itkTransformFileReader.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkRigid3DTransform.h>
#include <itkTransformFactory.h>
#include "itkTranslationTransform.h"
#include "itkVersorRigid3DTransform.h"

#include "GetPot.h"


namespace itk
{

  ApplyTransformToTensorImageCommand::ApplyTransformToTensorImageCommand()
  {
    m_ShortDescription = "Apply a rigid transformation to a tensor image without resampling";
    m_LongDescription = m_ShortDescription;
    m_LongDescription = "\n\nUsage:\n";
    m_LongDescription +="-i [input image]\n";
    m_LongDescription +="-t [transform]\n";
    m_LongDescription +="-o [Output file name]\n";
  }

  ApplyTransformToTensorImageCommand::~ApplyTransformToTensorImageCommand()
  {}

  int ApplyTransformToTensorImageCommand::Execute (int narg, const char* arg[])
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
  
    typedef double ScalarType;  
    typedef itk::Image<ScalarType, 3>                      ImageType;
    typedef itk::TensorImageIO<ScalarType, 3, 3>           IOType;
    typedef IOType::TensorImageType                        TensorImageType;
    typedef TensorImageType::PixelType                     TensorType;
    typedef TensorImageType::DirectionType DirectionType;
  
    typedef itk::MatrixOffsetTransformBase<double, 3, 3> LinearTransformType;
    typedef itk::VersorRigid3DTransform<double> RigidTransformType;
    typedef itk::TranslationTransform<double, 3> TranslationTransformType;
    typedef itk::Rigid3DTransform<double> Rigid3DTransformType;
    typedef itk::Transform<double, 3, 3> TransformType;
    typedef itk::TransformFileReader TransformReaderType;
    typedef TensorImageType::PointType PointType;
  
    itk::TransformFactory< LinearTransformType >::RegisterTransform ();
    itk::TransformFactory< RigidTransformType >::RegisterTransform ();
    itk::TransformFactory< TranslationTransformType >::RegisterTransform ();
    itk::TransformFactory< Rigid3DTransformType >::RegisterTransform ();
  
    std::cout<<"reading input file "<<fileIn<<std::endl;  
    IOType::Pointer myReader = IOType::New();
    myReader->SetFileName (fileIn);
    try
    {
      myReader->Read();
    } catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return EXIT_FAILURE;
    }  
    TensorImageType::Pointer myTensorImage = myReader->GetOutput();
    std::cout<<"done."<<std::endl;
  
    TransformReaderType::Pointer reader = TransformReaderType::New();
    reader->SetFileName ( fileTr );
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return EXIT_FAILURE;
    }
  
    LinearTransformType::Pointer transform1 = dynamic_cast<LinearTransformType*>( reader->GetTransformList()->front().GetPointer() );
    TranslationTransformType::Pointer transform2 = dynamic_cast<TranslationTransformType*>( reader->GetTransformList()->front().GetPointer() );
    Rigid3DTransformType::Pointer transform3 = dynamic_cast<Rigid3DTransformType*>( reader->GetTransformList()->front().GetPointer() );

    if (!transform1 && !transform2 && !transform3)
    {
      std::cerr << "The transformation written in "<<fileTr<<" does not derive from a MatrixOffsetTransformBase, "
		<<"which is mandatory in this executable to be able to only change the header of the file"<<std::endl;
      return EXIT_FAILURE;
    }
  
    PointType origin = myTensorImage->GetOrigin();
    PointType neworigin;
    if (transform1)
      neworigin = transform1->TransformPoint (origin);
    else if (transform2)
      neworigin = transform2->TransformPoint (origin);
    else if (transform3)
      neworigin = transform3->TransformPoint (origin);
    
    DirectionType direction = myTensorImage->GetDirection();
    DirectionType newdirection;
    // if (transform1)
    //   newdirection = direction * transform1->GetMatrix();
    // else if (transform2)
    //   newdirection = direction;
    // else if (transform3)
    //   newdirection = direction * transform1->GetMatrix();
    if (transform1)
      newdirection = transform1->GetMatrix() * direction;
    else if (transform2)
      newdirection = direction;
    else if (transform3)
      newdirection = transform3->GetMatrix() * direction;
  
    myTensorImage->SetOrigin(neworigin);
    myTensorImage->SetDirection(newdirection);
  
    std::cout<<"writing output file "<<fileOut<<std::endl;  
    IOType::Pointer myWriter = IOType::New();
    myWriter->SetFileName (fileOut);
    myWriter->SetInput (myTensorImage);
  
    try
    {
      myWriter->Write();
    } catch (itk::ExceptionObject &e)
    {
      std::cerr << e;
      return -1;
    }

  
    std::cout << " Done." << std::endl;  
  
    return EXIT_SUCCESS;

  }
}
