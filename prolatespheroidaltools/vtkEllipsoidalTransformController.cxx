#include <vtkEllipsoidalTransformController.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkLandmarkWidget.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPropPicker.h>
#include <vtkCellPicker.h>
#include <vtkPicker.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkInformationVector.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRendererCollection.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include <itkImageFileWriter.h>

vtkCxxRevisionMacro(vtkEllipsoidalTransformController, "$Revision: 1315 $");
vtkStandardNewMacro(vtkEllipsoidalTransformController);

//----------------------------------------------------------------------------
vtkEllipsoidalTransformControllerCommand::vtkEllipsoidalTransformControllerCommand()
{
  this->Controller = NULL;
}
//----------------------------------------------------------------------------
vtkEllipsoidalTransformControllerCommand::~vtkEllipsoidalTransformControllerCommand()
{
  if (this->Controller)
    this->Controller->UnRegister(this);
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformControllerCommand::SetController (vtkEllipsoidalTransformController* arg)
{
  if (arg == this->Controller)
    return;
  if (this->Controller)
    this->Controller->UnRegister(this);
  this->Controller = arg;
  if (this->Controller)
    this->Controller->Register(this);
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformControllerCommand::Execute ( vtkObject *caller, unsigned long event, void *vtkNotUsed(callData))
{

  if (!this->Controller)
    return;
  if (!this->Controller->GetInteractorCollection())
    return;
  
  if ( event == vtkCommand::EndInteractionEvent )
  {
    /***********  Landmark has finished moving  ***********/
    this->Controller->VerifyOrthogonality();
    this->Controller->RefreshConstraints();
    this->Controller->Update();
    
    this->Controller->GetInteractorCollection()->InitTraversal();
    vtkRenderWindowInteractor* interactor = vtkRenderWindowInteractor::SafeDownCast (this->Controller->GetInteractorCollection()->GetNextItemAsObject());
    while(interactor)
    {
      interactor->Render();
      interactor = vtkRenderWindowInteractor::SafeDownCast (this->Controller->GetInteractorCollection()->GetNextItemAsObject());
    }
  }
}

//----------------------------------------------------------------------------
vtkEllipsoidalTransformController::vtkEllipsoidalTransformController()
{
  this->SetNumberOfInputPorts(0);
  m_Input                    = NULL;
  this->InteractorCollection = NULL;
  this->Enabled              = 0;
  m_DomainFilter             = DomainFilterType::New();
  m_DomainEllipsoidFilter    = DomainEllipsoidFilterType::New();
  m_Transform                = TransformType::New();

  PointType p1, p2, p3, p4;

  p1[0] = 0; p1[1] = 0; p1[2] = 0;
  p2[0] = 0; p2[1] = 0; p2[2] = 3;
  p3[0] = 2; p3[1] = 0; p3[2] = 0;
  p4[0] = 0; p4[1] = 1; p4[2] = 0;
  m_Transform->DefineEllipsoid (p1,p2,p3,p4);
  
  this->WallMaxAngle = 95;
  this->WallThickness = 13;
  
  this->LandmarkCollection = vtkCollection::New();
  this->Command            = vtkEllipsoidalTransformControllerCommand::New();
  this->TotalLandmarkCollection = vtkCollection::New();
  this->Command->SetController (this);

  m_DomainEllipsoidFilter->SetInput (m_DomainFilter->GetOutput());
  m_DomainEllipsoidFilter->SetTransform (m_Transform);
}

//----------------------------------------------------------------------------
vtkEllipsoidalTransformController::~vtkEllipsoidalTransformController()
{
  this->LandmarkCollection->Delete();
  this->Command->Delete();
  if (this->InteractorCollection)
    this->InteractorCollection->UnRegister(this);
  this->TotalLandmarkCollection->Delete();
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::SetInteractorCollection (vtkCollection* arg)
{
  if (arg == this->InteractorCollection)
    return;
  if (this->InteractorCollection)
    this->InteractorCollection->UnRegister(this);
  this->InteractorCollection = arg;
  if (this->InteractorCollection)
    this->InteractorCollection->Register(this);
  this->SetEnabled (this->GetEnabled());
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::SetInput(ImageType::Pointer input)
{
  m_Input = input;
  m_DomainFilter->SetInput (input);
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::SetTransform (TransformType::Pointer arg)
{
  m_Transform = arg;
  
  if (this->GetInteractorCollection())
  {
    this->RemoveAllConstraints();
    
    TransformType::ParametersType points = arg->GetParameters();
    
    PointType p1, p2, p3, p4;
    p1[0] = points[0]; p1[1] = points[1]; p1[2] = points[2];
    p2[0] = points[3]; p2[1] = points[4]; p2[2] = points[5];
    p3[0] = points[6]; p3[1] = points[7]; p3[2] = points[8];
    p4[0] = points[9]; p4[1] = points[10]; p4[2] = points[11];
    
    this->AddConstraint (p1.GetDataPointer(), 0);
    this->AddConstraint (p2.GetDataPointer(), 1);
    this->AddConstraint (p3.GetDataPointer(), 2);
    this->AddConstraint (p4.GetDataPointer(), 3);
    
    this->RefreshConstraints();
  }

  m_DomainEllipsoidFilter->SetTransform (arg);

  this->Modified();
}

void vtkEllipsoidalTransformController::SetTransform (vtkPointSet* arg)
{
  if (!arg->GetPoints()->GetNumberOfPoints() || (arg->GetNumberOfPoints() != 4) )
  {
    vtkErrorMacro (<<"vtkEllipsoidalTransformController::SetConstraints() : argument not compatible");
    return;
  }
  
  TransformType::Pointer transform = TransformType::New();
  PointType p1 = arg->GetPoint (0);
  PointType p2 = arg->GetPoint (1);
  PointType p3 = arg->GetPoint (2);
  PointType p4 = arg->GetPoint (3);
  transform->DefineEllipsoid (p1,p2,p3,p4);
  this->SetTransform (transform);
}


//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::GuessTransform (vtkPointSet* arg)
{
  double percentile = 0.85;
  
  int numberofpoints = arg->GetNumberOfPoints();
  if (!numberofpoints)
  {
    std::cerr<<"vtkEllipsoidalTransformController::GuessTransform: there is no point in the argument, skipping"<<std::endl;
    TransformType::Pointer transform = TransformType::New();  
    this->SetTransform (transform);
    return;
  }
  
  /// compute the ventricle's center of mass
  itk::Point<double,3> point;
  itk::Point<double,3> mean; mean[0] = mean[1] = mean[2] = 0.0;
  for (int i=0; i<numberofpoints; i++)
  {
    arg->GetPoint (i, point.GetDataPointer());
    mean += point.GetVectorFromOrigin();
  }
  for (int i=0; i<3; i++) mean[i] /= (double)numberofpoints;
  
  /// compute the covariance matrix of all ventricle's points
  vnl_matrix_fixed<double,3,3> matrix (0.0);
  for (int i=0; i<numberofpoints; i++)
  {
    arg->GetPoint (i, point.GetDataPointer());
    for (unsigned int j=0; j<3; j++)
      for (unsigned int k=0; k<3; k++)
	matrix[j][k] += (point[j] - mean[j])*(point[k] - mean[k]);
  }
  matrix /= (double)numberofpoints;
  
  itk::Vector<double,3> longaxis;
  itk::Vector<double,3> shortaxis1, shortaxis2;
  typedef vnl_symmetric_eigensystem< double >  SymEigenSystemType;
  SymEigenSystemType eig (matrix);
  for(unsigned int i=0;i<3;i++)
    longaxis[i] = eig.V(i,2);
  for(unsigned int i=0;i<3;i++)
    shortaxis1[i] = eig.V(i,1);
  for(unsigned int i=0;i<3;i++)
    shortaxis2[i] = eig.V(i,2);
  
  longaxis.Normalize();
  shortaxis1.Normalize();
  shortaxis2.Normalize();
  
  /// compute coordinates of apex and base
  double projectionrange[2] = {VTK_DOUBLE_MAX, VTK_DOUBLE_MIN};
  itk::Point<double,3> apex, base;
  apex[0] = apex[1] = apex[2] = 0.0;
  base[0] = base[1] = base[2] = 0.0;
  
  itk::Vector<double,3> vec;
  for (int i=0; i<numberofpoints; i++)
  {
    arg->GetPoint (i, point.GetDataPointer());
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
    arg->GetPoint (i, point.GetDataPointer());
    vec = point - mean;
    double projection = vec * longaxis;
    if (projection < middle)
      continue;
    radius += std::sqrt ( vec.GetSquaredNorm() - projection * projection );
    midnumberofpoints ++;
  }
  radius /= (double)midnumberofpoints;
  itk::Point<double,3> septum1 = base + radius * shortaxis1;
  itk::Point<double,3> septum2 = base + radius * shortaxis2;
  
  /// Fill the Prolate Spheroidal Transform with information
  TransformType::ParametersType parameters;
  parameters.SetSize (12);
  parameters[0] = base[0];   parameters[1] = base[1];   parameters[2] = base[2];
  parameters[3] = apex[0];   parameters[4] = apex[1];   parameters[5] = apex[2]; 
  parameters[6] = septum1[0]; parameters[7] = septum1[1]; parameters[8] = septum1[2];
  parameters[9] = septum2[0]; parameters[10] = septum2[1]; parameters[11] = septum2[2];
  TransformType::Pointer transform = TransformType::New();
  transform->SetParameters (parameters);
  
  this->SetTransform (transform);
}

//----------------------------------------------------------------------------
vtkLandmarkWidget* vtkEllipsoidalTransformController::AddConstraint (double* pos, int type)
{
  if (!this->InteractorCollection)
    return NULL;
  this->GetInteractorCollection()->InitTraversal();
  vtkRenderWindowInteractor* item = vtkRenderWindowInteractor::SafeDownCast (this->GetInteractorCollection()->GetNextItemAsObject());
  vtkLandmarkWidget* initial_landmark = NULL;
  while(item)
  {
    vtkLandmarkWidget* l = vtkLandmarkWidget::New();
    l->ScaleOff();
    l->SetCenter (pos);
    l->SetRadius (3.0);
    l->SetValue (type);
    if (type == 0) l->GetSphereProperty()->SetColor (1,0,0);
    if (type == 1) l->GetSphereProperty()->SetColor (1,1,0);
    if (type == 2) l->GetSphereProperty()->SetColor (0,1,0);
    if (type == 3) l->GetSphereProperty()->SetColor (0,0,1);
    l->SetRepresentationToSurface();
    l->SetInteractor (item);
    l->SetEnabled (true);
    this->GetTotalLandmarkCollection()->AddItem (l);
    if (!initial_landmark)
      initial_landmark = l;
    l->Delete();
    item = vtkRenderWindowInteractor::SafeDownCast (this->GetInteractorCollection()->GetNextItemAsObject());
  }
  this->LandmarkCollection->AddItem(initial_landmark);
  this->LinkInteractions();
  return initial_landmark;
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::RemoveAllConstraints (void)
{
  this->LandmarkCollection->InitTraversal();
  vtkLandmarkWidget* landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetNextItemAsObject());
  while(landmark)
  {
    landmark->Off();
    landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetNextItemAsObject());
  }
  this->LandmarkCollection->RemoveAllItems();
}
  
//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::RefreshConstraints (void)
{
  vtkLandmarkWidget* landmark = NULL;
  landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetItemAsObject(0));
  PointType p1 (landmark->GetCenter());
  landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetItemAsObject(1));
  PointType p2 (landmark->GetCenter());
  landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetItemAsObject(2));
  PointType p3 (landmark->GetCenter());
  landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetItemAsObject(3));
  PointType p4 (landmark->GetCenter());
  
  m_Transform->DefineEllipsoid (p1,p2,p3,p4);

  m_DomainEllipsoidFilter->SetTransform (m_Transform);

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::GetLandmarkSet(vtkPolyData* arg)
{
  vtkPoints* pts = vtkPoints::New();
  vtkFloatArray* values = vtkFloatArray::New();
  this->LandmarkCollection->InitTraversal();
  vtkLandmarkWidget* landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetNextItemAsObject());
  while(landmark)
  {
    pts->InsertNextPoint (landmark->GetCenter());
    values->InsertNextValue (landmark->GetValue());
    landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetNextItemAsObject());
  }
  arg->SetPoints (pts);
  arg->GetPointData()->SetScalars (values);
  pts->Delete();
  values->Delete();
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::SetEnabled (unsigned int arg)
{
  this->Enabled = arg;

  this->LandmarkCollection->InitTraversal();
  vtkLandmarkWidget* landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetNextItemAsObject());
  while(landmark)
  {
    landmark->SetEnabled(arg);
    landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetNextItemAsObject());
  }
}

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::LinkInteractions ( void)
{ 
  vtkLandmarkWidget* l1 = NULL;
  vtkLandmarkWidget* l2 = NULL;
  unsigned int numberofinteractors = this->GetInteractorCollection()->GetNumberOfItems();
  vtkCollection* collection        = this->GetTotalLandmarkCollection();
  unsigned int numberoflandmarks   = collection->GetNumberOfItems() / numberofinteractors;
  for (unsigned int n = 0; n < numberoflandmarks; n++)
  {
    for (unsigned int i=0; i<numberofinteractors; i++)
    {
      int id1 = numberofinteractors * n + i;
      l1 = vtkLandmarkWidget::SafeDownCast (collection->GetItemAsObject(id1));
      if (!l1->HasObserver(vtkCommand::EndInteractionEvent, this->Command))
	l1->AddObserver(vtkCommand::EndInteractionEvent, this->Command, -1);
      
      for (unsigned int j=0; j<numberofinteractors; j++)
      {
	int id2 = numberofinteractors * n + j;
	l2 = vtkLandmarkWidget::SafeDownCast (collection->GetItemAsObject(id2));
	
	if (l1 != l2)
	{
	  if (!l1->HasObserver(vtkCommand::InteractionEvent, l2->GetCommand()))
	    l1->AddObserver(vtkCommand::InteractionEvent, l2->GetCommand());
	  if (!l1->HasObserver(vtkCommand::EnableEvent, l2->GetCommand()))
	    l1->AddObserver(vtkCommand::EnableEvent,      l2->GetCommand());
	  if (!l1->HasObserver(vtkCommand::DisableEvent, l2->GetCommand()))
	    l1->AddObserver(vtkCommand::DisableEvent,     l2->GetCommand());
	  if (!l2->HasObserver(vtkCommand::InteractionEvent, l1->GetCommand()))
	    l2->AddObserver(vtkCommand::InteractionEvent, l1->GetCommand());
	  if (!l2->HasObserver(vtkCommand::EnableEvent, l1->GetCommand()))
	    l2->AddObserver(vtkCommand::EnableEvent,      l1->GetCommand());
	  if (!l2->HasObserver(vtkCommand::DisableEvent, l1->GetCommand()))
	    l2->AddObserver(vtkCommand::DisableEvent,     l1->GetCommand());
	}  
      }    
    }
  }
}


//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::CreateGrid (TransformType::Pointer transform,
						    double lambda1, double lambda2,
						    double mu1, double mu2,
						    double nu1, double nu2,
						    unsigned int throughwalldivisions,
						    unsigned int longdivisions,
						    unsigned int circumdivisions,
						    vtkPolyData* grid,
						    bool stay_in_prolate,
						    bool imitatecircular)
{
  unsigned int N[3], N1[3], N2[3];
  N[0] = throughwalldivisions + 1;
  N[1] = longdivisions + 1;
  N[2] = circumdivisions + 1;
  
  unsigned int start_long_axis = 0;
  unsigned int stop_long_axis = N[1];
  unsigned int start_circum = 0;
  unsigned int stop_circum = N[2];
  
  N1[0] = 0;
  N1[1] = start_long_axis;
  N1[2] = start_circum;
  
  N2[0] = N[0];
  N2[1] = stop_long_axis;
  N2[2] = stop_circum;
  
  double bounds[3][2];
  bounds[0][0] = lambda1; bounds[0][1] = lambda2;
  bounds[1][0] = mu1;     bounds[1][1] = mu2;
  bounds[2][0] = nu1;     bounds[2][1] = nu2;
  
  vtkPoints* outputpoints = vtkPoints::New();
  
  grid->SetPoints (outputpoints);
  grid->Allocate();

  double zeta1factor = 40;
  double zeta2factor = 4;

  for (int octant = 1; octant < 8; octant+=2)
    for (unsigned int i=N1[0]; i<N2[0]; i++)
      for (unsigned int j=N1[1]; j<N2[1]; j++)
	for (unsigned int k=N1[2]; k<N2[2]; k++)
	{
	  vtkIdType cell[2];
	  
	  TransformType::InputPointType point1, point2, pt;
	  point1[0] = bounds[0][0] + (double)i * ( (bounds[0][1] - bounds[0][0]) / (double)(N[0] - 1) );
	  point1[1] = bounds[1][0] + (double)j * ( (bounds[1][1] - bounds[1][0]) / (double)(N[1] - 1) );
	  point1[2] = bounds[2][0] + (double)k * ( (bounds[2][1] - bounds[2][0]) / (double)(N[2] - 1) );
	  point1[3] = octant;
	  
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
	    point2[3] = octant;
	    
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
	    point2[3] = octant;
	    
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
	    point2[3] = octant;
	    
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

//----------------------------------------------------------------------------
void vtkEllipsoidalTransformController::VerifyOrthogonality (void)
{
  vtkLandmarkWidget* landmark = NULL, *l3 = NULL, *l4 = NULL;
  landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetItemAsObject(0));
  PointType p1 (landmark->GetCenter());
  landmark = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetItemAsObject(1));
  PointType p2 (landmark->GetCenter());
  l3 = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetItemAsObject(2));
  PointType p3 (l3->GetCenter());
  l4 = vtkLandmarkWidget::SafeDownCast (this->LandmarkCollection->GetItemAsObject(3));
  PointType p4 (l4->GetCenter());
  
  VectorType v1 = p2 - p1;
  VectorType v2 = p3 - p1;
  VectorType v3 = p4 - p1;
  if (v2.GetNorm() >= v1.GetNorm())
  {
    v2.Normalize();
    p3 = p1 + (v1.GetNorm() -  vcl_numeric_limits<double>::epsilon()) * v2;
  }
  if (v3.GetNorm() >= v2.GetNorm())
  {
    v3.Normalize();
    p4 = p1 + (v2.GetNorm() -  vcl_numeric_limits<double>::epsilon()) * v3;
  }
  
  v1.Normalize();
  v2.Normalize();
  v3 = CrossProduct (v1,v2);
  v2 = CrossProduct (v3,v1);
  double shortaxis1 = p1.EuclideanDistanceTo (p3);
  double shortaxis2 = p1.EuclideanDistanceTo (p4);
  p3 = p1 + shortaxis1 * v2;
  p4 = p1 + shortaxis2 * v3;
  l3->SetCenter (p3.GetDataPointer());
  l4->SetCenter (p4.GetDataPointer());
}

//----------------------------------------------------------------------------
int vtkEllipsoidalTransformController::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  unsigned int circumdivisions = 20;
  unsigned int longdivisions = 15;
  unsigned int throughwalldivisions = 3;
  //  double maxangle = this->WallMaxAngle;
  double thickness = this->WallThickness;
  
  double a2 = m_Transform->GetLambda1() * m_Transform->GetLambda1();
  double b2 = m_Transform->GetLambda2() * m_Transform->GetLambda2();
  double c2 = m_Transform->GetLambda3() * m_Transform->GetLambda3();
  
  double lambda1 = 0.0 - 20 * thickness;
  double lambda2 = 0.0 + 20 * thickness;    
  double mu1 = c2;
  double mu2 = b2;
  double nu1 = b2;
  double nu2 = a2;
  
  if ( mu2 <= mu1)
  {
    vtkGenericWarningMacro (<<"ellipsoid transform ill-posed with Lambda3 >= Lambda2 = "<<m_Transform->GetLambda2()<<"\n");
    mu2 = mu1 + 1.0;
  }
  
  if ( nu2 <= nu1)
  {
    vtkGenericWarningMacro (<<"ellipsoid transform ill-posed with Lambda2 >= Lambda1 = "<<m_Transform->GetLambda1()<<"\n");
    nu2 = nu1 + 1.0;
  }
  
  TransformType::Pointer inversetransform = TransformType::New();
  m_Transform->GetInverse (inversetransform);
  
  PointType basalpoint = m_Transform->GetCenter();
  VectorType longaxis  = m_Transform->GetAxis1();
  
  m_DomainFilter->SetCroppingPlane (basalpoint, -longaxis);
  m_DomainEllipsoidFilter->SetCroppingPlane (basalpoint, -longaxis);
  m_DomainEllipsoidFilter->SetWallThickness (this->WallThickness);
  m_DomainEllipsoidFilter->SetWallMaxAngle (this->WallMaxAngle);
  
  m_DomainFilter->Update();
  // m_DomainEllipsoidFilter->Update();
  
  this->CreateGrid (inversetransform,
		    lambda1, lambda2,
		    mu1, mu2,
		    nu1, nu2,
		    throughwalldivisions,
		    longdivisions,
		    circumdivisions,
		    output);
  return 1;
}
