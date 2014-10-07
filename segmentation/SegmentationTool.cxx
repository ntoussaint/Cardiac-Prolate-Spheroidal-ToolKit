/*=========================================================================

Module:    $Id: SegmentationTest.cxx 917 2008-08-27 10:37:34Z ntoussaint $
Language:  C++
Author:    $Author: ntoussaint $
Date:      $Date: 2008-08-27 12:37:34 +0200 (Wed, 27 Aug 2008) $
Version:   $Revision: 917 $

=========================================================================*/

#include <vtksys/SystemTools.hxx>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageActor.h>
#include <vtkCamera.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataReader.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkDataSetWriter.h>
#include <vtkImageData.h>
#include <vtkCommand.h>
#include <vtkProperty.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkExtractImageFilter.h>

#include <vtkImageViewCollection.h>
#include <vtkImageView2D.h>
#include <vtkImageView3D.h>
#include <vtkImageViewCornerAnnotation.h>
#include <vtkImageMapToColors.h>
#include <vtkInteractorStyleImageView2D.h>
#include <vtkImageViewCollection.h>

#include <vtkLandmarkSegmentationController.h>

#include "GetPot.h"



typedef vtkLandmarkSegmentationController::ScalarType ScalarType;
typedef vtkLandmarkSegmentationController::ImageType  ImageType;
typedef itk::Image<ScalarType, 4>                     Image4DType;
typedef itk::ImageFileReader<ImageType>               ReaderType;
typedef itk::ImageFileReader<Image4DType>             Reader4DType;
typedef itk::ImageToVTKImageFilter<ImageType>         ConverterType;
typedef itk::ImageToVTKImageFilter<Image4DType>       Converter4DType;
typedef ImageType::DirectionType                      DirectionType;
typedef ImageType::PointType                          PointType;


void PrintHelp(const std::string exec)
{
  std::cout << std::endl << exec.c_str() << " Usage: " << std::endl << std::endl;
  std::cout << "-i [input image file]" << std::endl;
  std::cout << "-a [extra input image file for visualization, repeat if necessary]" << std::endl;
  std::cout << "-l [optional initial landmark-set]" << std::endl;
  std::cout << "-e [optional extra mesh file]" << std::endl;
  std::cout << "-o [output file base (use 'epi' or 'endo']" << std::endl << std::endl;

  std::ostringstream help;
  help << "Please take notice of the interactive usage : \n"
       << "This tool is dedicated to the interactive creation of a segmentation mesh\n"
       << "from a set of points interactively placed in space.\n"
       << "You can add/remove/move points using the following on any open window :\n"
       << "\t -Shift+MiddleClick    : adds a point that will lie on the target surface.\n"
       << "\t -Shift+LeftClick      : adds a point that will lie in the interior of the target surface.\n"
       << "\t -Shift+RightClick     : adds a point that will lie in the exterior of the target surface.\n"
       << "\t -Ctrl+Shift+LeftClick : removes the targetted point.\n"
       << "\t -LeftClick & move     : drags   the targetted point.\n\n"
       << "*** IMPORTANT *** : Once the segmentation is satisfactory, Save the results by pressing 's' in the 3D view.\n";
       
  std::cout << std::endl << help.str().c_str() << std::endl;
  
  std::exit(EXIT_SUCCESS);
}

// -----------------------------------------------------------------------------
class vtkMyCommand : public vtkCommand
{
public:
  static vtkMyCommand *New() {return new vtkMyCommand;};
  void SetSegmenter (vtkLandmarkSegmentationController* arg)
  { this->Segmenter = arg; };
  void SetViewer (vtkImageView3D* arg)
  { this->Viewer = arg; };
  void SetFileOut (const char* arg)
  { this->FileOut = arg; }
  
  virtual void Execute(vtkObject *caller,unsigned long event, void *vtkNotUsed(callData))
  {
    if (event == vtkCommand::CharEvent)
    {
      vtkRenderWindowInteractor* rwi = vtkRenderWindowInteractor::SafeDownCast(caller);
      if (rwi->GetKeyCode() == 'v')
      {
        this->Viewer->SetRenderingMode (!this->Viewer->GetRenderingMode());
	this->Viewer->Render();
      }
      
      if (rwi->GetKeyCode() == 's')
      {
	if (this->Segmenter)
	{
	  std::ostringstream meshfile;
	  meshfile << this->FileOut << "-mesh.vtk";
	  std::ostringstream landmarkfile;
	  landmarkfile << this->FileOut << "-landmarks.vtk";
	  std::ostringstream functionfile;
	  functionfile << this->FileOut << "-function.mha";
	  
	  std::cout<<"writing "<<meshfile.str().c_str()<<" ... "<<std::flush;
	  vtkDataSetWriter* writer = vtkDataSetWriter::New();
	  writer->SetInputConnection (this->Segmenter->GetOutputPort());
	  writer->SetFileName (meshfile.str().c_str());	
	  writer->Update();
	  std::cout<<"done."<<std::endl;
	  std::cout<<"writing "<<landmarkfile.str().c_str()<<" ... "<<std::flush;
	  vtkPolyData* lms = vtkPolyData::New();
	  this->Segmenter->GetLandmarkSet(lms);
	  writer->SetInputData (lms);
	  writer->SetFileName (landmarkfile.str().c_str());
	  writer->Update();
	  writer->Delete();
	  lms->Delete();
	  std::cout<<"done."<<std::endl;
	  std::cout<<"writing "<<functionfile.str().c_str()<<" ... "<<std::flush;
	  typedef itk::ImageFileWriter<vtkLandmarkSegmentationController::ImageType> WriterType;
	  WriterType::Pointer imagewriter = WriterType::New();
	  imagewriter->SetInput (this->Segmenter->GetImplicitImage());
	  imagewriter->SetFileName (functionfile.str().c_str());
	  imagewriter->Update();
	  std::cout<<"done."<<std::endl;
	}
      }
    }    
  }
 protected:
  vtkMyCommand(){};
  ~vtkMyCommand(){};
 private:
  vtkLandmarkSegmentationController* Segmenter;
  vtkImageView3D* Viewer;
  const char* FileOut;
};


ImageType::Pointer ExtractFirstItem (Image4DType::Pointer input)
{
  typedef Image4DType::RegionType Region4dType;
  typedef Image4DType::SpacingType Spacing4Dtype;
  typedef itk::ImageRegionIterator<Image4DType> Iterator4DType;
  typedef Iterator4DType::IndexType Index4DType;
  typedef ImageType::DirectionType Direction3Dtype;
  typedef Image4DType::DirectionType Direction4Dtype;
  typedef itk::ExtractImageFilter<Image4DType, ImageType> ExtractImageType;
  Image4DType::RegionType region = input->GetLargestPossibleRegion();
  Image4DType::SizeType   size   = region.GetSize();
  Image4DType::IndexType  index;
  index[0] = 0; index[1] = 0; index[2] = 0; index[3] = 0;
  size[3] = 0;
  region.SetSize (size);
  region.SetIndex (index);
  ExtractImageType::Pointer extractor = ExtractImageType::New();
  extractor->SetExtractionRegion (region);
  extractor->SetInput (input);
  extractor->SetDirectionCollapseToSubmatrix();
  
  try
  {
    extractor->Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr<<e<<std::endl;
    std::exit (EXIT_FAILURE);
  }
  
  return extractor->GetOutput();
}


// -----------------------------------------------------------------------------
int main (int argc, char* argv[])
{

  GetPot   cl(argc, argv); // argument parser
  if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    PrintHelp(cl[0]);

  const bool IsInputPresent     = cl.search(2,"-i","-I");
  const bool IsOutputPresent    = cl.search(2,"-o","-O");
  const bool IsLandmarksPresent = cl.search(2,"-l","-L");
  const bool IsExtraPresent     = cl.search(2,"-e","-E");
  
  if (!IsInputPresent || !IsOutputPresent)
  {
    std::cerr << "Error: Input and/or output not set." << std::endl;
    std::exit (EXIT_FAILURE);
  }
  
  const std::string fileIn        = cl.follow("NoFile",2,"-i","-I");
  const std::string fileExtra     = cl.follow("NoFile",2,"-e","-E");
  const std::string fileOut       = cl.follow("NoFile",2,"-o","-O");
  const std::string fileLandmarks = cl.follow("NoFile",2,"-l","-L");

  std::vector<std::string> imagefiles;  
  imagefiles.push_back (std::string (fileIn));
  
  // read all arguments that start with '-a'
  cl.init_multiple_occurrence();
  std::string file = cl.follow((const char*)"", "-a");
  while(file.size() != 0) {
    imagefiles.push_back(std::string(file));
    file = cl.follow((const char*)"", "-a");
  }
  
  int position[2] = {0, 0};
  vtkImageViewCollection*      pool = vtkImageViewCollection::New();
  vtkImageView3D*            view3d = vtkImageView3D::New();
  vtkRenderWindowInteractor* iren3d = vtkRenderWindowInteractor::New();
  vtkRenderWindow*           rwin3d = vtkRenderWindow::New();
  vtkRenderer*                ren3d = vtkRenderer::New();
  
  position[0] = 500; position[1] = 90;
  iren3d->SetRenderWindow(rwin3d);
  rwin3d->AddRenderer (ren3d);
  view3d->SetRenderWindow(rwin3d);
  view3d->SetRenderer(ren3d);
  view3d->SetPosition (position);
  double color[3]={0.9,0.9,0.9};
  view3d->SetBackground (color);
  
  vtkMyCommand* command = vtkMyCommand::New();
  view3d->GetInteractor()->AddObserver(vtkCommand::CharEvent, command);
  pool->AddItem (view3d);
  command->SetViewer (view3d);
  
  ImageType::Pointer firstimage = NULL;
  vtkImageData* firstinput = NULL;
  vtkMatrix4x4* firstmatrix = NULL;

  vtkProperty* p1 = vtkProperty::New();
  p1->LightingOff();
  p1->SetColor (1,0,0);
  vtkProperty* p2 = vtkProperty::New();
  p2->LightingOff();
  p2->SetColor (0,1,0);
  vtkProperty* p3 = vtkProperty::New();
  p3->LightingOff();
  p3->SetColor (0,0,1);
  
  std::vector<ConverterType::Pointer> listofconverters;
  
  for (unsigned int N=0; N<imagefiles.size(); N++)
  {
    std::cout<<"reading "<<imagefiles[N].c_str()<<" ..."<<std::endl;
    Reader4DType::Pointer reader = Reader4DType::New();
    reader->SetFileName (imagefiles[N].c_str());
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr<<e<<std::endl;
      std::exit (EXIT_FAILURE);
    }
    
    ImageType::Pointer image = ExtractFirstItem (reader->GetOutput());
    PointType origin = image->GetOrigin();
    ImageType::DirectionType direction = image->GetDirection();
    
    ConverterType::Pointer converter = ConverterType::New();
    converter->SetInput (image);
    try
    {
      converter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr<<e<<std::endl;
      std::exit (EXIT_FAILURE);
    }
    vtkImageData* input = converter->GetOutput();
    
    listofconverters.push_back (converter);
    
    vtkMatrix4x4* matrix = vtkMatrix4x4::New();
    matrix->Identity();
    
    for (unsigned int i=0; i<3; i++)
      for (unsigned int j=0; j<3; j++)
	matrix->SetElement (i,j,direction[i][j]);
    PointType correctedorigin;
    matrix->MultiplyPoint (origin.GetDataPointer(), correctedorigin.GetDataPointer());
    for (int i=0; i<3; i++)
      matrix->SetElement (i, 3, origin[i] - correctedorigin[i]);    
    
    vtkImageView2D* view = vtkImageView2D::New();
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    vtkRenderWindow* rwin = vtkRenderWindow::New();
    vtkRenderer* ren = vtkRenderer::New();
    position[0] = 130; position[1] = 90;
    iren->SetRenderWindow(rwin);
    rwin->AddRenderer (ren);
    view->SetRenderWindow(rwin);
    view->SetRenderer(ren);
    view->SetPosition (position);
    view->SetInput (input);
    view->SetOrientationMatrix(matrix);
    view->GetCornerAnnotation()->SetText (1, imagefiles[N].c_str());

    pool->AddItem (view);
    view3d->AddExtraPlane (view->GetImageActor());
    view->Delete();
    iren->Delete();
    ren->Delete();
    rwin->Delete();

    if (N == 0)
    {
      firstimage  = image;
      firstinput  = input;
      firstmatrix = matrix;
      view->SetViewOrientation (vtkImageView2D::VIEW_ORIENTATION_AXIAL);
      
      vtkImageView2D* view2 = vtkImageView2D::New();
      vtkRenderWindowInteractor* iren2 = vtkRenderWindowInteractor::New();
      vtkRenderWindow* rwin2 = vtkRenderWindow::New();
      vtkRenderer* ren2 = vtkRenderer::New();
      position[0] = 500; position[1] = 480;
      iren2->SetRenderWindow(rwin2);
      rwin2->AddRenderer (ren2);
      view2->SetRenderWindow(rwin2);
      view2->SetRenderer(ren2);
      view2->SetPosition (position);
      view2->SetInput (input);
      view2->SetOrientationMatrix(matrix);
      view2->SetViewOrientation (vtkImageView2D::VIEW_ORIENTATION_SAGITTAL);
      view2->GetCornerAnnotation()->SetText (1, fileIn.c_str());
      
      pool->AddItem (view2);
      int* sl_r = view2->GetSliceRange();
      view2->SetSlice ( (int)( (sl_r[1] - sl_r[1]) / 2.0) );
      view3d->AddExtraPlane (view2->GetImageActor());
  
      view2->Delete();
      iren2->Delete();
      rwin2->Delete();
      ren2->Delete();
      
      vtkImageView2D* view3 = vtkImageView2D::New();
      vtkRenderWindowInteractor* iren3 = vtkRenderWindowInteractor::New();
      vtkRenderWindow* rwin3 = vtkRenderWindow::New();
      vtkRenderer* ren3 = vtkRenderer::New();
      position[0] = 130; position[1] = 480;
      iren3->SetRenderWindow(rwin3);
      rwin3->AddRenderer (ren3);
      view3->SetRenderWindow(rwin3);
      view3->SetRenderer(ren3);
      view3->SetPosition (position);
      view3->SetInput (input);
      view3->SetOrientationMatrix(matrix);
      view3->SetViewOrientation (vtkImageView2D::VIEW_ORIENTATION_CORONAL);
      view3->GetCornerAnnotation()->SetText (1, fileIn.c_str());

      pool->AddItem (view3);    
      view3d->AddExtraPlane (view3->GetImageActor());
      view3->Delete();
      iren3->Delete();
      rwin3->Delete();
      ren3->Delete();
    }
  }
  
  view3d->SetInput (firstinput);
  view3d->SetOrientationMatrix(firstmatrix);
  
  pool->SyncReset();
  pool->SyncSetShowAnnotations (1);
  pool->SyncSetShowImageAxis (1);
  pool->SyncSetShowImageAxis (1);
  pool->SyncSetShowScalarBar (0);
  pool->SyncSetInterpolate (0);
  pool->SyncSetAnnotationStyle (vtkImageView2D::AnnotationStyle2);
  pool->SyncSetViewConvention (vtkImageView2D::VIEW_CONVENTION_RADIOLOGICAL);
  pool->LinkColorWindowLevelOn ();
  pool->LinkSliceMoveOn();
  pool->SyncSetWheelInteractionStyle(vtkInteractorStyleImageView2D::InteractionTypeSlice);
  
  if (IsExtraPresent)
  {
    std::cout<<"reading "<<fileExtra<<" ..."<<std::endl;
    vtkDataSetReader* r = vtkDataSetReader::New();
    r->SetFileName (fileExtra.c_str());
    r->Update();
    if (vtkPointSet::SafeDownCast (r->GetOutput()))
      pool->SyncAddDataSet (vtkPointSet::SafeDownCast (r->GetOutput()));
    r->Delete();
  }
  
  // --------------------------------------------------------
  // ------------ Part concerning segmentation --------------
  // --------------------------------------------------------
  
  vtkCollection* interactorcollection = vtkCollection::New();
  pool->InitTraversal();
  vtkImageView* item = pool->GetNextItem();
  while(item) { interactorcollection->AddItem (item->GetInteractor()); item = pool->GetNextItem(); }
  vtkLandmarkSegmentationController* segmentation = vtkLandmarkSegmentationController::New();
  segmentation->SetInteractorCollection (interactorcollection);

  segmentation->SetInput (firstimage);
  segmentation->EnabledOn();
  
  if (IsLandmarksPresent)
  {
    std::cout<<"reading "<<fileLandmarks<<std::endl;
    vtkPolyDataReader* r = vtkPolyDataReader::New();
    r->SetFileName (fileLandmarks.c_str());
    r->Update();
    segmentation->SetConstraints (r->GetOutput());
    segmentation->Update();
    r->Delete();
  }
  
  pool->SyncAddDataSet (segmentation->GetOutput());

  // --------------------------------------------------------
  // ---------- End Part concerning segmentation ------------
  // --------------------------------------------------------

  command->SetSegmenter (segmentation);
  command->SetFileOut (fileOut.c_str());
  
  int size[2]={370, 370};
  pool->SyncSetSize (size);
  pool->SyncRender();  
  pool->SyncStart();
  view3d->Delete();
  iren3d->Delete();
  rwin3d->Delete();
  ren3d->Delete();
  pool->Delete();
  command->Delete();
  interactorcollection->Delete();
  segmentation->Delete();
  
  return EXIT_SUCCESS;
}
