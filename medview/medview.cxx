/*=========================================================================

Module:    $Id: SegmentationTest.cxx 917 2008-08-27 10:37:34Z ntoussaint $
Language:  C++
Author:    $Author: ntoussaint $
Date:      $Date: 2008-08-27 12:37:34 +0200 (Wed, 27 Aug 2008) $
Version:   $Revision: 917 $

=========================================================================*/

#include <vtksys/SystemTools.hxx>
#include <vtkRenderer.h>
#include <vtkDataArrayCollection.h>
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

#include <vtkDataManager.h>
#include <vtkMetaImageData.h>
#include <vtkMetaDataSetSequence.h>
#include <vtkMetaSurfaceMesh.h>
#include <vtkMetaVolumeMesh.h>

#include "vtkContourWidget.h"
#include "vtkOrientedGlyphContourRepresentation.h"

#include "GetPot.h"

void PrintHelp(const std::string exec)
{
  std::cout << std::endl << exec.c_str() << " Usage: " << std::endl << std::endl;
  std::cout <<  exec.c_str() << " input-image-1 [input-image-2] [...] [-m <multiple-view>" << std::endl;
  
  std::exit(EXIT_SUCCESS);
}


// -----------------------------------------------------------------------------
class vtkMyCommand : public vtkCommand
{
  
public:
  
  static vtkMyCommand *New() {return new vtkMyCommand;};

  void SetViewer (vtkImageView3D* viewer)
  { this->Viewer = viewer; };
  void SetCollection (vtkImageViewCollection* arg)
  { this->Collection = arg; };
  void AddSequence (vtkMetaDataSetSequence* sequence)
  { this->Sequences.push_back (sequence); };
  
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

      if (rwi->GetKeyCode() == 'y')
      {
      	if (!this->Sequences.size())
      	  return;
      	
	int id = this->Sequences[0]->GetCurrentId();
	if (id >= this->Sequences[0]->GetNumberOfMetaDataSets() - 1)
	  id = 0;
	else
	  id = id+1;
	for (unsigned int i=0; i<this->Sequences.size(); i++)
	  this->Sequences[i]->UpdateToIndex (id);
	this->Collection->SyncRender();
	
      }
      if (rwi->GetKeyCode() == 'z')
      {
	this->Collection->SyncSetViewOrientation (orientation++);
	if (orientation > 2) orientation = 0;
	this->Collection->SyncRender();
      }
    }
  }
  

 protected:
  vtkMyCommand()
  {
    toggle = 1;
    orientation = 0;
    this->Viewer = NULL;
    this->Collection = NULL;
  };
  ~vtkMyCommand(){};

 private:
  vtkImageView3D* Viewer;
  vtkDataManager* Manager;
  vtkImageViewCollection* Collection;
  
  std::vector<vtkMetaDataSetSequence*> Sequences;
  bool toggle;
  unsigned int orientation;
};



// -----------------------------------------------------------------------------
int main (int argc, char* argv[])
{

  GetPot   cl(argc, argv); // argument parser
  if( cl.size() == 1 || cl.search(2, "--help", "-h") )
    PrintHelp(cl[0]);

  const bool MultiView          = cl.follow((const bool)(0), "-m");
  
  if (argc < 2)
  {
    std::cerr << "Error: Input not set." << std::endl;
    std::exit (EXIT_FAILURE);
  }
  
  std::vector<std::string> imagefiles;
  // read all arguments
  for (int i=1; i<argc; i++)
  {
    std::string file = argv[i];
    if (file.size() > 3)
      imagefiles.push_back(file);
  }
  
  int position[2] = {0, 0};
  vtkImageViewCollection*      pool = vtkImageViewCollection::New();
  vtkImageView3D*            view3d = vtkImageView3D::New();
  vtkRenderWindowInteractor* iren3d = vtkRenderWindowInteractor::New();
  vtkRenderWindow*           rwin3d = vtkRenderWindow::New();
  vtkRenderer*                ren3d = vtkRenderer::New();

  vtkDataManager*           manager = vtkDataManager::New();
  
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
  command->SetCollection (pool);
  
  vtkImageData* firstinput = NULL;
  vtkMatrix4x4* firstmatrix = NULL;
  vtkImageView2D* firstview = NULL;
  
  for (unsigned int N=0; N<imagefiles.size(); N++)
  {
    std::cout<<"reading "<<imagefiles[N].c_str()<<" ..."<<std::endl;

    vtkMetaDataSet* metadataset = NULL;
    
    try
    {
      metadataset = manager->ReadFile(imagefiles[N].c_str());
    } catch (...)
    {
      // not a metadataset...
    }

    if ( metadataset->GetType() == vtkMetaDataSet::VTK_META_IMAGE_DATA )
    {
      vtkMetaImageData* metaimage = vtkMetaImageData::SafeDownCast(metadataset);
      vtkImageData* image = NULL;
      
      if (!metaimage)
      {
	metaimage = vtkMetaImageData::SafeDownCast(vtkMetaDataSetSequence::SafeDownCast(metadataset)->GetMetaDataSet(0));
	image = vtkImageData::SafeDownCast(vtkMetaDataSetSequence::SafeDownCast(metadataset)->GetDataSet());
      }
      else
	image = metaimage->GetImageData();
      
      vtkMatrix4x4* matrix = metaimage->GetOrientationMatrix();

      vtkImageData* input = image;
      
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

      vtkOrientedGlyphContourRepresentation *contourRep = vtkOrientedGlyphContourRepresentation::New();
      vtkContourWidget *contourWidget = vtkContourWidget::New();
      contourWidget->SetInteractor(iren);
      contourWidget->SetRepresentation(contourRep);
      contourWidget->SetKeyPressActivation (1);
      contourWidget->SetKeyPressActivationValue ('j');
      contourRep->Delete();
      
      view->Delete();
      iren->Delete();
      ren->Delete();
      rwin->Delete();
      
      if (N == 0)
      {
	firstinput  = input;
	firstmatrix = matrix;
	firstview   = view;
      }
      
      if ( (N == 0) && (MultiView) )
      {
	
	view->SetSliceOrientation (vtkImageView2D::SLICE_ORIENTATION_XZ);
	
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
	view2->SetSliceOrientation (vtkImageView2D::SLICE_ORIENTATION_YZ);
	view2->GetCornerAnnotation()->SetText (1, imagefiles[N].c_str());
	
	pool->AddItem (view2);
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
	view3->SetSliceOrientation (vtkImageView2D::SLICE_ORIENTATION_XY);
	view3->GetCornerAnnotation()->SetText (1, imagefiles[N].c_str());
	
	pool->AddItem (view3);
	view3d->AddExtraPlane (view3->GetImageActor());
	
	view3->Delete();
	iren3->Delete();
	rwin3->Delete();
	ren3->Delete();
      }
    }
    else
    {
      vtkProperty* prop = vtkProperty::SafeDownCast( metadataset->GetProperty() );
      prop->SetColor (0.5,0.5,0.5);
      prop->SetLineWidth (2);
      // if (vtkPointSet::SafeDownCast (metadataset->GetDataSet())->GetNumberOfPoints() > 10000)
      // 	prop->SetOpacity (0.9);
      pool->SyncAddDataSet( vtkPointSet::SafeDownCast (metadataset->GetDataSet()));
      metadataset->AddActor (vtkActor::SafeDownCast (view3d->FindDataSetActor (metadataset->GetDataSet())));
      std::cout<<*(metadataset->GetDataSet())<<std::endl;
      
      vtkDataArrayCollection* arrays = vtkDataArrayCollection::New();
      metadataset->GetColorArrayCollection(arrays);
      arrays->InitTraversal();
      vtkDataArray* item = arrays->GetNextItem();
      metadataset->ColorByArray (item);
      while(item)
      {
	std::cout<<"mesh has array called : "<<item->GetName()<<std::endl;
	item = arrays->GetNextItem();
      }
      
      
    }
    
    if (vtkMetaDataSetSequence::SafeDownCast(metadataset))
      command->AddSequence (vtkMetaDataSetSequence::SafeDownCast (metadataset));
  }    
  
  if (firstinput && firstmatrix)
  {
    view3d->SetInput (firstinput);
    view3d->SetOrientationMatrix(firstmatrix);
  }
  
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
  manager->Delete();
  
  return EXIT_SUCCESS;
}
