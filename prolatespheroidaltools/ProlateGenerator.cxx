#include <vtksys/SystemTools.hxx>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleImage.h>
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
#include <itkTransformFileWriter.h>
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

#include <vtkImageViewCollection.h>
#include <vtkImageView2D.h>
#include <vtkImageView3D.h>
#include <vtkImageViewCornerAnnotation.h>
#include <vtkInteractorStyleImageView2D.h>
#include <vtkImageViewCollection.h>

#include <vtkProlateSpheroidalTransformController.h>
#include <itkProlateSpheroidalTransform.h>

#include "GetPot.h"

void PrintHelp(const char* exec)
{
  std::cout << std::endl << exec << " Usage: " << std::endl << std::endl;
  std::cout << "-i  [input image file]" << std::endl;
  std::cout << "-l1 [landmark-set 1 (default : endo-landmarks.vtk)]" << std::endl;
  std::cout << "-l2 [landmark-set 2 (default : epi-landmarks.vtk)]" << std::endl;
  std::cout << "-p  [optional initial prolate-transform file]" << std::endl;
  std::cout << "-wt [optional wall thickness in mm (default 13)]" << std::endl;
  std::cout << "-wa [optional wall basal max angle in deg (default 95)]" << std::endl;
  std::cout << "-o  [output file base (default: '')]" << std::endl << std::endl;
  
  std::exit(EXIT_SUCCESS);
}

// -----------------------------------------------------------------------------
class vtkMyCommand : public vtkCommand
{
public:

  typedef vtkLandmarkSegmentationController::ScalarType ScalarType;
  typedef vtkLandmarkSegmentationController::ImageType  ImageType;
  typedef vtkLandmarkSegmentationController::BooleanImageType BooleanImageType;
  typedef itk::ImageFileReader<ImageType>               ReaderType;
  typedef itk::ImageToVTKImageFilter<ImageType>         ConverterType;
  typedef ImageType::DirectionType                      DirectionType;
  typedef ImageType::PointType                          PointType;
  typedef ImageType::SpacingType                        VectorType;
  typedef itk::ProlateSpheroidalTransform<ScalarType>        TransformType;
  
  static vtkMyCommand *New() {return new vtkMyCommand;};
  void SetViewer (vtkImageView3D* arg)
  { this->Viewer = arg; };
  void SetFileOut (const char* arg)
  { this->FileOut = arg; this->FileOutSet = true;}
  void SetController (vtkProlateSpheroidalTransformController* arg)
  { this->Controller = arg; }
  
  virtual void Execute(vtkObject *caller,unsigned long event, void *vtkNotUsed(callData))
  {
    if (event == vtkCommand::CharEvent)
    {
      vtkRenderWindowInteractor* rwi = vtkRenderWindowInteractor::SafeDownCast(caller);
      if (rwi->GetKeyCode() == 's')
      {
	std::ostringstream meshfile1;
	if (this->FileOutSet)
	  meshfile1 << this->FileOut;
	meshfile1 << "domain.vtk";
	std::ostringstream imagefile1;
	if (this->FileOutSet)
	  imagefile1 << this->FileOut;
	imagefile1 << "domain.mha";
	std::ostringstream transformfile;
	if (this->FileOutSet)
	  transformfile << this->FileOut;
	transformfile << "prolatetransform.tr";
	std::ostringstream imagefile2;
	if (this->FileOutSet)
	  imagefile2 << this->FileOut;
	imagefile2 << "domain-ellipsoid.mha";
	std::ostringstream meshfile2;
	if (this->FileOutSet)
	  meshfile2 << this->FileOut;
	meshfile2 << "domain-ellipsoid.vtk";

	std::cout<<"writing "<<imagefile1.str().c_str()<<" ... "<<std::flush;
	typedef itk::ImageFileWriter<BooleanImageType> WriterType;
	WriterType::Pointer imagewriter = WriterType::New();
	imagewriter->SetInput (this->Controller->GetDomain());
	imagewriter->SetFileName (imagefile1.str().c_str());
	imagewriter->Update();
	std::cout<<"done."<<std::endl;
	std::cout<<"writing "<<meshfile1.str().c_str()<<" ... "<<std::flush;
	vtkDataSetWriter* writer = vtkDataSetWriter::New();
	writer->SetInputData (this->Controller->GetDomainMesh());
	writer->SetFileName (meshfile1.str().c_str());
	writer->Update();
	std::cout<<"done."<<std::endl;
	std::cout<<"writing "<<imagefile2.str().c_str()<<" ... "<<std::flush;
	imagewriter->SetInput (this->Controller->GetDomainEllipsoid());
	imagewriter->SetFileName (imagefile2.str().c_str());
	imagewriter->Update();
	std::cout<<"done."<<std::endl;
	std::cout<<"writing "<<meshfile2.str().c_str()<<" ... "<<std::flush;
	writer->SetInputData (this->Controller->GetDomainEllipsoidMesh());
	writer->SetFileName (meshfile2.str().c_str());
	writer->Update();
	std::cout<<"done."<<std::endl;
	
	/// Write the transform
	std::cout<<"writing "<<transformfile.str().c_str()<<" ... "<<std::flush;
	itk::TransformFactory<TransformType>::RegisterTransform ();
	itk::TransformFileWriter::Pointer transformwriter = itk::TransformFileWriter::New();
	transformwriter->SetFileName(transformfile.str().c_str());
	transformwriter->SetInput (this->Controller->GetTransform());
	transformwriter->Update();
	std::cout<<"done."<<std::endl;

      }
    }
  }
 protected:
  vtkMyCommand()
  {};
  ~vtkMyCommand()
  {
    this->FileOut = NULL;
    this->FileOutSet = false;
  };
 private:
  vtkImageView3D* Viewer;
  vtkProlateSpheroidalTransformController* Controller;
  const char* FileOut;
  bool FileOutSet;
  
};

// -----------------------------------------------------------------------------
int main (int argc, char* argv[])
{

  GetPot   cl(argc, argv); // argument parser
  if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    PrintHelp(cl[0].c_str());

  const bool IsInputPresent     = cl.search(2,"-i","-I");
  const bool IsOutputPresent    = cl.search(2,"-o","-O");
  const bool IsTransformPresent = cl.search(2,"-p","-P");
  
  if (!IsInputPresent)
  {
    std::cerr << "Error: Input and/or output not set." << std::endl;
    std::exit (EXIT_FAILURE);
  }
  
  typedef vtkProlateSpheroidalTransformController::ScalarType ScalarType;
  typedef vtkProlateSpheroidalTransformController::ImageType  ImageType;
  typedef vtkProlateSpheroidalTransformController::TransformType  TransformType;
  typedef itk::ImageFileReader<ImageType>               ReaderType;
  typedef itk::ImageToVTKImageFilter<ImageType>         ConverterType;
  typedef ImageType::DirectionType                      DirectionType;
  typedef ImageType::PointType                          PointType;
  
  const std::string fileIn = cl.follow("NoFile",2,"-i","-I");
  const std::string fileOut        = cl.follow("",2,"-o","-O");
  const std::string fileLandmarks1 = cl.follow("endo-landmarks.vtk",2,"-l1","-L1");
  const std::string fileLandmarks2 = cl.follow("epi-landmarks.vtk",2,"-l2","-L2");
  const std::string transformfile = cl.follow("NoFile",2,"-p","-P");
  const double wallthickness = cl.follow(13, 2,"-wt","-WT");
  const double wallmaxangle = cl.follow(95, 2,"-wa","-WA");

  std::cout<<"reading "<<fileIn.c_str()<<std::endl;
  
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName (fileIn.c_str());
  reader->Update();
  ConverterType::Pointer converter = ConverterType::New();
  converter->SetInput (reader->GetOutput());
  converter->Update();
  vtkImageData* input = converter->GetOutput();
  
  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  matrix->Identity();
  DirectionType direction = reader->GetOutput()->GetDirection();
  PointType origin = reader->GetOutput()->GetOrigin();
  for (unsigned int i=0; i<3; i++)
    for (unsigned int j=0; j<3; j++)
      matrix->SetElement (i,j,direction[i][j]);
  PointType correctedorigin;
  matrix->MultiplyPoint (origin.GetDataPointer(), correctedorigin.GetDataPointer());
  for (int i=0; i<3; i++)
    matrix->SetElement (i, 3, origin[i]-correctedorigin[i]);

  int position[2] = {0, 0};
  vtkImageViewCollection*      pool = vtkImageViewCollection::New();
  vtkImageView3D*            view3d = vtkImageView3D::New();
  vtkRenderWindowInteractor* iren3d = vtkRenderWindowInteractor::New();
  vtkRenderWindow*           rwin3d = vtkRenderWindow::New();
  vtkRenderer*                ren3d = vtkRenderer::New();

  position[0] = 400; position[1] = 420;
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
  
  vtkImageView2D* view = vtkImageView2D::New();
  vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
  vtkRenderWindow* rwin = vtkRenderWindow::New();
  vtkRenderer* ren = vtkRenderer::New();
  position[0] = 0; position[1] = 0;
  iren->SetRenderWindow(rwin);
  rwin->AddRenderer (ren);
  view->SetRenderWindow(rwin);
  view->SetRenderer(ren);
  view->SetPosition (position);
  view->SetInput (input);
  view->SetOrientationMatrix(matrix);
  view->GetCornerAnnotation()->SetText (1, fileIn.c_str());
  view->SetViewOrientation (vtkImageView2D::VIEW_ORIENTATION_AXIAL);
  pool->AddItem (view);  
  // view3d->AddExtraPlane (view->GetImageActor());
  view->Delete();
  iren->Delete();
  ren->Delete();
  rwin->Delete();

  vtkImageView2D* view2 = vtkImageView2D::New();
  vtkRenderWindowInteractor* iren2 = vtkRenderWindowInteractor::New();
  vtkRenderWindow* rwin2 = vtkRenderWindow::New();
  vtkRenderer* ren2 = vtkRenderer::New();
  position[0] = 400; position[1] = 0;
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
  // view3d->AddExtraPlane (view2->GetImageActor());
  view2->Delete();
  iren2->Delete();
  rwin2->Delete();
  ren2->Delete();
  
  vtkImageView2D* view3 = vtkImageView2D::New();
  vtkRenderWindowInteractor* iren3 = vtkRenderWindowInteractor::New();
  vtkRenderWindow* rwin3 = vtkRenderWindow::New();
  vtkRenderer* ren3 = vtkRenderer::New();
  position[0] = 0; position[1] = 420;
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
  // view3d->AddExtraPlane (view3->GetImageActor());
  view3->Delete();
  iren3->Delete();
  rwin3->Delete();
  ren3->Delete();
  
  vtkImageView2D* firstview = vtkImageView2D::SafeDownCast (pool->GetItem (1));
  if (firstview)
  {
    view3d->SetInput (firstview->GetInput());
    view3d->SetOrientationMatrix(firstview->GetOrientationMatrix());
  }

  pool->SyncReset();
  pool->SyncSetShowAnnotations (1);
  pool->SyncSetShowImageAxis (1);
  pool->SyncSetShowScalarBar (0);
  pool->SyncSetAnnotationStyle (vtkImageView2D::AnnotationStyle2);
  pool->SyncSetViewConvention (vtkImageView2D::VIEW_CONVENTION_RADIOLOGICAL);
  pool->LinkColorWindowLevelOn ();
  pool->LinkSliceMoveOn();
  pool->SyncSetWheelInteractionStyle(vtkInteractorStyleImageView2D::InteractionTypeSlice);  

  // --------------------------------------------------------
  // ---- Part concerning Prolate Spheroidal Generation -----
  // --------------------------------------------------------
  
  vtkPolyDataReader* r1 = vtkPolyDataReader::New();
  r1->SetFileName (fileLandmarks1.c_str());
  r1->Update();

  vtkPolyDataReader* r2 = vtkPolyDataReader::New();
  r2->SetFileName (fileLandmarks2.c_str());
  r2->Update();

  vtkCollection* interactorcollection = vtkCollection::New();
  pool->InitTraversal();
  vtkImageView* item = pool->GetNextItem();
  while(item) { interactorcollection->AddItem (item->GetInteractor()); item = pool->GetNextItem(); }  
  vtkProlateSpheroidalTransformController* controller = vtkProlateSpheroidalTransformController::New();

  controller->SetInteractorCollection (interactorcollection);
  controller->SetInput (reader->GetOutput());
  controller->SetLandmarkSet1 (r1->GetOutput());
  controller->SetLandmarkSet2 (r2->GetOutput());
  controller->SetWallThickness (wallthickness);
  controller->SetWallMaxAngle (wallmaxangle);
  controller->EnabledOn();
  
  if (IsTransformPresent)
  {
    TransformType::Pointer transform = TransformType::New();
    itk::TransformFactory<TransformType>::RegisterTransform ();
    itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
    transformreader->SetFileName(transformfile.c_str());
    transformreader->Update();
    transform = dynamic_cast<TransformType*>( transformreader->GetTransformList()->front().GetPointer() );
    controller->SetTransform (transform);
  }
  else
  {
    controller->Update();
    controller->GuessTransform(controller->GetDomainMesh());
  }
  
  vtkProperty* prop1 = vtkProperty::New();
  prop1->SetColor (0.5,0.5,0.5);
  prop1->SetOpacity (0.3);
  vtkProperty* prop2 = vtkProperty::New();
  prop2->SetColor (0.6,0.0,0.0);
  prop2->SetOpacity (1);
  
  pool->SyncAddDataSet (controller->GetDomainMesh(), prop1);
  pool->SyncAddDataSet (controller->GetOutput(), prop2);  

  if (IsOutputPresent) command->SetFileOut (fileOut.c_str());
  command->SetController (controller);  
  
  // --------------------------------------------------------
  // -- End Part concerning Prolate Spheroidal Generation ---
  // --------------------------------------------------------
  

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
  prop1->Delete();
  prop2->Delete();

  return EXIT_SUCCESS;
}
