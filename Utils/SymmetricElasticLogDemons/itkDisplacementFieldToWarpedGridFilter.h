#ifndef _itkDisplacementFieldToWarpedGridFilter_h_
#define _itkDisplacementFieldToWarpedGridFilter_h_

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImageToMeshFilter.h"
#include "itkAutomaticTopologyMeshSource.h"

namespace itk
{

/**
 * @class DisplacementFieldToWarpedGridFilter
 *
 * @brief Create a grid mesh and warp it according to the input deformation
 *        field. Grid and warp can be 3D or 2D (with respect to a plane)
 */

/// 3D ONLY
template < class TImage, class TMesh >
class ITK_EXPORT DisplacementFieldToWarpedGridFilter :
    public ImageToMeshFilter< TImage, TMesh >
{

public:

  /// Usual typedefs
  typedef DisplacementFieldToWarpedGridFilter  Self;
  typedef ImageToMeshFilter< TImage, TMesh >   Superclass;
  typedef SmartPointer<Self>                   Pointer;
  typedef SmartPointer<const Self>             ConstPointer;

  /// Method for creation through the object factory
  itkNewMacro( Self );

  /// Run-time type information (and related methods)
  itkTypeMacro( DisplacementFieldToWarpedGridFilter, ImageToMeshFilter );

  /// Type of available grids */
  enum GridType {
     GRID_3D = 0,
     GRID_2D_X,
     GRID_2D_Y,
     GRID_2D_Z
  };

protected:
  typedef itk::Vector< double, 3 >     VectorType;

  typedef TImage                       ImageType;
  typedef typename TImage::Pointer     ImagePointer;
  typedef typename TImage::PixelType   PixelType;
  typedef typename TImage::PointType   ImagePointType;

  typedef TMesh                        MeshType;
  typedef typename TMesh::PointType    MeshPointType;

  typedef typename itk::AutomaticTopologyMeshSource< MeshType >  MeshSourceType;
  typedef typename MeshSourceType::Pointer                       MeshSourcePointer;

public:

  /// Set the type of grid
  virtual void SetGridType( unsigned int gtype )
  { m_GridType = GridType(gtype); }

  /// Set/Get grid resolution
  itkSetMacro( GridResolution, unsigned int );
  itkGetMacro( GridResolution, unsigned int );

/*   /// Set/Get inverse input deformation field */
/*   itkSetMacro( InverseInputField, bool ); */
/*   itkGetMacro( InverseInputField, bool ); */
/*   itkBooleanMacro( InverseInputField ); */

protected:

  /// Default constructor
  DisplacementFieldToWarpedGridFilter();

  /// Default destructor
  ~DisplacementFieldToWarpedGridFilter();

  /// Generate the data: elastic filtering of the input image
  virtual void GenerateData( void );

  /// Generate stacks of plane along a given direction
  virtual void GeneratePlanes(
      MeshSourcePointer meshSource, MeshPointType origin,
      VectorType step1, VectorType step2, VectorType position );

/*   /// Invert input field */
/*   virtual ImagePointer InvertInputField( void ); */

private:

  DisplacementFieldToWarpedGridFilter(Self&);  // intentionally not implemented
  void operator=(const Self&);                 // intentionally not implemented

  /// Inverse input deformation field
  bool m_InverseInputField;

  /// Grid type (3D, 2D_X, 2D_Y, 2D_Z)
  GridType m_GridType;

  /// Grid resolution
  unsigned int m_GridResolution;
};


}   // end namespace itk




#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDisplacementFieldToWarpedGridFilter.txx"
#endif

#endif
