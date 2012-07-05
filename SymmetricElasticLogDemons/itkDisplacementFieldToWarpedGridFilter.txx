#ifndef _itkDisplacementFieldToWarpedGridFilter_txx
#define _itkDisplacementFieldToWarpedGridFilter_txx

#include "itkDisplacementFieldToWarpedGridFilter.h"

#include "itkImageToMeshFilter.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkAutomaticTopologyMeshSource.h"
#include "itkMesh.h"
#include "itkTriangleCell.h"
#include "itkWarpMeshFilter.h"
#include "itkVector.h"

namespace itk
{


// -----------------------------------------------------------------------------
// Constructor
template < class TImage, class TMesh >
DisplacementFieldToWarpedGridFilter< TImage, TMesh >
::DisplacementFieldToWarpedGridFilter()
{
  m_GridResolution = 50;
  m_GridType = GRID_3D;
//   m_InverseInputField = true;
}



// -----------------------------------------------------------------------------
// Destructor
template < class TImage, class TMesh >
DisplacementFieldToWarpedGridFilter< TImage, TMesh >
::~DisplacementFieldToWarpedGridFilter()
{
}



// -----------------------------------------------------------------------------
// Generate planes

template < class TImage, class TMesh >
void
DisplacementFieldToWarpedGridFilter< TImage, TMesh >
::GeneratePlanes( MeshSourcePointer meshSource, MeshPointType origin,
    VectorType step1, VectorType step2, VectorType position )
{
  if ( !meshSource ) // Should never happen
    return;

  MeshPointType p1, p2, p3;
  for ( float plane = 0; plane <= m_GridResolution; ++plane )
    for ( float i = 1.0; i <= m_GridResolution; ++i )
      for ( float j = 1.0; j <= m_GridResolution; ++j )
      {
        p1 = origin + step1 * (i-1)/m_GridResolution + step2 * (j-1)/m_GridResolution + position * plane/m_GridResolution;
        p2 = origin + step1 *  i/m_GridResolution    + step2 * (j-1)/m_GridResolution + position * plane/m_GridResolution;
        p3 = origin + step1 *  i/m_GridResolution    + step2 *  j/m_GridResolution    + position * plane/m_GridResolution;
        meshSource->AddTriangle(
            meshSource->AddPoint( p1 ),
            meshSource->AddPoint( p2 ),
            meshSource->AddPoint( p3 )
        );

        p1 = origin + step1 * (i-1)/m_GridResolution + step2 * (j-1)/m_GridResolution + position * plane/m_GridResolution;
        p2 = origin + step1 * (i-1)/m_GridResolution + step2 *  j/m_GridResolution    + position * plane/m_GridResolution;
        p3 = origin + step1 *  i/m_GridResolution    + step2 *  j/m_GridResolution    + position * plane/m_GridResolution;
        meshSource->AddTriangle(
            meshSource->AddPoint( p1 ),
            meshSource->AddPoint( p2 ),
            meshSource->AddPoint( p3 )
        );
      }

}


// // -----------------------------------------------------------------------------
// // Invert input deformation field using def inverse

// template < class TImage, class TMesh >
// typename DisplacementFieldToWarpedGridFilter< TImage, TMesh>::ImagePointer
// DisplacementFieldToWarpedGridFilter< TImage, TMesh >
// ::InvertInputField( void )
// {
//   // Do not work for image dimensions different from 3
//   if ( ImageType::ImageDimension != 3 )
//     exit( EXIT_FAILURE );

//   typename ImageType::IndexType idx;

//   typename ImageType::RegionType::SizeType size =
//     this->GetInput(0)->GetLargestPossibleRegion().GetSize();

//   // We split the input ITK vector image into N inrimage scalar images
//   // Pasha definverse only works with float images
//   yav::Inrimage* fieldinr =
//     new yav::Inrimage( size[0], size[1], size[2],
//         yav::Inrimage::WT_FLOAT_VECTOR, 3, VM_NON_INTERLACED );
//   fieldinr->allocate();
//   float**** fieldarray = (float****) (fieldinr->getArray());

//   typedef typename itk::ImageRegionConstIteratorWithIndex< ImageType > ConstIteratorType;
//   ConstIteratorType itt( this->GetInput(0), this->GetInput(0)->GetLargestPossibleRegion() );
//   for ( itt.GoToBegin(); !itt.IsAtEnd(); ++itt )
//   {
//     PixelType t = itt.Get();

//     idx = itt.GetIndex();
//     fieldarray[0][idx[2]][idx[1]][idx[0]] = static_cast<float>(t[0]);
//     fieldarray[1][idx[2]][idx[1]][idx[0]] = static_cast<float>(t[1]);
//     fieldarray[2][idx[2]][idx[1]][idx[0]] = static_cast<float>(t[2]);
//   }

//   yav::Inrimage** splittedField = fieldinr->splitVectors();

//   // Invert deformation field
//   yav::Transfo* t = new yav::Transfo(
//       splittedField[0], splittedField[1], splittedField[2],
//       splittedField[0]->getX(), splittedField[0]->getY(), splittedField[0]->getZ());
//   t->inverse();

//   float*** x = (float***) t->x->getArray();
//   float*** y = (float***) t->y->getArray();
//   float*** z = (float***) t->z->getArray();

//   // Convert back the inverted deformation field to an ITK image
//   ImagePointer invertedField = ImageType::New();
//   invertedField->CopyInformation( this->GetInput(0) );
//   invertedField->SetRegions( this->GetInput(0)->GetLargestPossibleRegion() );
//   invertedField->Allocate();
//   typedef typename itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
//   IteratorType it( invertedField, invertedField->GetLargestPossibleRegion() );
//   PixelType p;
//   for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
//   {
//     idx = it.GetIndex();
//     p[0] = x[idx[2]][idx[1]][idx[0]];
//     p[1] = y[idx[2]][idx[1]][idx[0]];
//     p[2] = z[idx[2]][idx[1]][idx[0]];
//     it.Set( p );
//   }

//   return invertedField;
// }


// -----------------------------------------------------------------------------
// GenerateData

template < class TImage, class TMesh >
void
DisplacementFieldToWarpedGridFilter< TImage, TMesh >
::GenerateData()
{
  // Get image bounding box in real coordinate system
  typename ImageType::RegionType::SizeType size =
    this->GetInput(0)->GetLargestPossibleRegion().GetSize();

  typename ImageType::IndexType idx;
  ImagePointType origin, point;

  VectorType dx, dy, dz;

  idx[0] = 0; idx[1] = 0; idx[2] = 0;
  this->GetInput(0)->TransformIndexToPhysicalPoint( idx, origin );

  idx[0] = size[0]-1; idx[1] = 0; idx[2] = 0;
  this->GetInput(0)->TransformIndexToPhysicalPoint( idx, point );
  dx = point - origin;

  idx[0] = 0; idx[1] = size[1]-1; idx[2] = 0;
  this->GetInput(0)->TransformIndexToPhysicalPoint( idx, point );
  dy = point - origin;

  idx[0] = 0; idx[1] = 0; idx[2] = size[2]-1;
  this->GetInput(0)->TransformIndexToPhysicalPoint( idx, point );
  dz = point - origin;

  // Create the grid
  // I have to do this weird copy to make the code compile on Mac OS X Tiger.
  // Apparently, the compiler cannot cast ImagePointType to MeshPointType...
  MeshPointType meshOrigin;
  for ( unsigned int i = 0; i < ImageType::ImageDimension; ++i )
    meshOrigin[i] = origin[i];

  MeshSourcePointer meshSource = MeshSourceType::New();
  if ( m_GridType == GRID_3D || m_GridType != GRID_2D_X )
    GeneratePlanes( meshSource, meshOrigin, dy, dz, dx );
  if ( m_GridType == GRID_3D || m_GridType != GRID_2D_Y )
    GeneratePlanes( meshSource, meshOrigin, dx, dz, dy );
  if ( m_GridType == GRID_3D || m_GridType != GRID_2D_Z )
    GeneratePlanes( meshSource, meshOrigin, dx, dy, dz );

  // Remove the third component of the deformation field if 2D deformation is
  // wanted
  typedef typename itk::ImageDuplicator< ImageType > DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( this->GetInput(0) );
  duplicator->Update();

  ImagePointer field = duplicator->GetOutput();

//   // Invert input deformation field if required
//   if ( m_InverseInputField )
//     field = this->InvertInputField();

//   if ( m_GridType != GRID_3D )
//   {
//     unsigned int componentToRemove = 0;
//     switch( m_GridType )
//     {
//     case GRID_2D_X: componentToRemove = 0; break;
//     case GRID_2D_Y: componentToRemove = 1; break;
//     case GRID_2D_Z: componentToRemove = 2; break;
//     default: break; // Should never happen
//     }

//     typedef typename itk::ImageRegionIterator< ImageType > IteratorType;
//     typename ImageType::PixelType pixel;
//     IteratorType it( field, field->GetLargestPossibleRegion() );
//     for ( it.GoToBegin(); ! it.IsAtEnd(); ++it )
//     {
//       pixel = it.Get();
//       pixel[componentToRemove] = 0;
//       it.Set( pixel );
//     }
//   }

  // Warp the grid
  typedef typename itk::WarpMeshFilter< MeshType, MeshType, ImageType > WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput( meshSource->GetOutput() );
  warper->SetDeformationField( field );
  warper->Update();

  // Save the result
  this->GraftOutput( warper->GetOutput() );
  this->GetOutput()->Modified();
}



} // end namespace itk

#endif

