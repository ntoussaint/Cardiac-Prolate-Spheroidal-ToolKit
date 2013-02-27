function(CreateCpstkCommand command_name keyword description)
  set(COMMAND_NAME  "itk${command_name}Command")
  set(FACTORY_NAME  "itk${command_name}CommandFactory")
  
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${COMMAND_NAME}.h ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.h ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.cxx
    COMMAND ${CMAKE_COMMAND} 
    -D SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cmake
    -D MY_COMMAND_NAME=${command_name}
    -D MY_DESCRIPTION=${description} 
    -D MY_KEYWORD=${keyword}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CreateCPSTKCommand.cmake
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/itkCpstkCommandTemplate.h.in
	    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/itkCpstkCommandFactoryTemplate.h.in
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/itkCpstkCommandFactoryTemplate.cxx.in
	    )
	  
  set(CPSTK_COMMAND_SRCS
    ${CPSTK_COMMAND_SRCS}
    ${COMMAND_NAME}.cxx
    ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.cxx
    PARENT_SCOPE
    )

endfunction()

CreateCpstkCommand(RotateProlateSpheroid
  "rotate" "Rotate a Prolate Spheroid according to a vtk file describing the antero-posterior line"
  )
CreateCpstkCommand(ExtractProlateInformation
  "extract" "Extract meaningful information in prolate spheroidal coordinates"
  )
CreateCpstkCommand(TensorImageToMesh
  "image2mesh" "Converts a tensor image into a tensor mesh structure"
  )
CreateCpstkCommand(TensorMeshToImage
  "mesh2image" "Converts a tensor mesh structure into a tensor image"
  )
CreateCpstkCommand(ExtrapolateTensorField
  "extrapolate" "Tensor Field Dense Approximation from Sparse Data"
  )
CreateCpstkCommand(CreateSyntheticCardiacTensorMap
  "synthesize" "Create a synthetic tensor map out of some cardiac LV fibre structure information"
  )
CreateCpstkCommand(FindClosestProlateTransform
  "find" "Find a prolate spheroid approximating an LV segmentation mesh"
  )
CreateCpstkCommand(CreateProlateDomain
  "domain" "Create an ellipsoidal domain image out of a prolate sph. transform"
  )
CreateCpstkCommand(TensorMeshStatistics
  "stats" "Compute statistics in prolate sph. coordinates over a set of tensors"
  )
CreateCpstkCommand(CumulateProlateSpheroidalDataSets
  "cumulate" "Cumulate several datasets in one, and normalizing 1st component"
  )
CreateCpstkCommand(LimitImageToAHAZone
  "aha-image" "Crop an image to an AHA zone in the Prolate Spheroidal sense"
  )
CreateCpstkCommand(LimitTensorsToAHAZone
  "aha-tensors" "Crop a tensor image or mesh to an AHA zone in the Prolate Spheroidal sense"
  )
CreateCpstkCommand(ApplyTransformToImage
  "apply-image" "Apply a rigid transformation to an image without resampling"
  )
CreateCpstkCommand(ApplyTransformToTensorImage
  "apply-tensors" "Apply a rigid transformation to tensor image without resampling"
  )
CreateCpstkCommand(ApplyTransformToMesh
  "apply-mesh" "Apply an affine transformation to a vtk mesh"
  )
CreateCpstkCommand(ResampleImage3
  "res-image" "Resample a 3D image to a reference field-of-view"
  )
CreateCpstkCommand(ResampleTensorImage3
  "res-tensors" "Resample a 3D tensor image to a reference field-of-view"
  )
CreateCpstkCommand(ColorFibresWithProlateAngle
  "colorify" "Colorify a fibre field with Prolate angular information"
  )
CreateCpstkCommand(OptimizeExtrapolationKernels
  "optimize" "Find optimal extrapolation Prol. Sph. kernel sizes to minimize discrepancies in a LS sense"
  )
CreateCpstkCommand(ExtractKernelsEnveloppe
  "enveloppe" "Extract the kernel distance enveloppes in n AHA zones"
  )
CreateCpstkCommand(TensorMeshGradient
  "gradient" "Compute the gradient of a tensor field in Prolate Spheroidal sense"
  )
CreateCpstkCommand(TensorMeshStructure
  "structure" "Compute the structure tensor field of a tensor field in Prolate Spheroidal sense"
  )
CreateCpstkCommand(TensorMeshCovariance
  "covariance" "Compute the covariance tensor field of a tensor field in Prolate Spheroidal sense"
  )
CreateCpstkCommand(ForwardTransformMesh
  "forward" "Transform a tensor image or mesh in Prolate Spheroidal Coordinates"
  )
CreateCpstkCommand(ForwardTransformImage
  "image2pss" "Transform a scalar image in a Prolate Spheroidal Coordinates box"
  )
CreateCpstkCommand(BackwardTransformImage
  "pss2image" "Transform Prolate Spheroidal Coordinates box into a scalar image"
  )
CreateCpstkCommand(BackwardTransformMesh
  "backward" "Transform a tensor image or mesh from Prolate Spheroidal Coordinates"
  )
CreateCpstkCommand(TensorsToVTK
  "itk2vtk" "Convert a tensor image or a list of tensor images into a vtkUnstructuredGrid structure"
  )
CreateCpstkCommand(TensorsToVTK2
  "itk2vtk2" "Convert a list of tensor images into a list vtkUnstructuredGrid structures"
  )
CreateCpstkCommand(CreateProlateCoordinateImage
  "coordinates" "Extract Prolate Coordinates from input and output it in an 3-compartment image"
  )
CreateCpstkCommand(ReorderDWIs
  "reorder" "Reorders DWIs in a coherent manner starting with the b-0"
  )
CreateCpstkCommand(CropImage
  "crop" "crop a 3D or 4D image to a certain diameter around the central long axis"
  )
CreateCpstkCommand(PaintTensorImageWithAngle
  "paint" "paint"
  )

add_library(CPSTKCommands ${CPSTK_BUILD_TYPE}
  ${CPSTK_COMMAND_SRCS}
  )

target_link_libraries(CPSTKCommands
  ${ITK_LIBRARIES}
  CPSTKProgramFactory
  CPSTKNumerics
  vtkCommon
  vtkIO
  ITKTensor
  )

if(NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
  install(TARGETS CPSTKCommands
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif(NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)