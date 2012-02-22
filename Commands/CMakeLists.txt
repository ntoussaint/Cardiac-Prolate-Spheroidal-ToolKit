function(CreatePtkCommand command_name keyword description)
  set(COMMAND_NAME  "itk${command_name}Command")
  set(FACTORY_NAME  "itk${command_name}CommandFactory")
  
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${COMMAND_NAME}.h ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.h ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.cxx
    COMMAND ${CMAKE_COMMAND} 
    -D SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cmake
    -D MY_COMMAND_NAME=${command_name}
    -D MY_DESCRIPTION=${description} 
    -D MY_KEYWORD=${keyword}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CreatePTKCommand.cmake
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/itkPtkCommandTemplate.h.in
	    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/itkPtkCommandFactoryTemplate.h.in
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/itkPtkCommandFactoryTemplate.cxx.in
	    )
	  
  set(PTK_COMMAND_SRCS
    ${PTK_COMMAND_SRCS}
    ${COMMAND_NAME}.cxx
    ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.cxx
    PARENT_SCOPE
    )
  
endfunction()



CreatePtkCommand(TensorsToVTK
  "itk2vtk" "Convert a tensor image or a list of tensor images into a vtkUnstructuredGrid structure"
  )
CreatePtkCommand(RotateProlateSpheroid
  "rotate" "Rotate a Prolate Spheroid according to a vtk file describing the antero-posterior line"
  )
CreatePtkCommand(ExtractProlateInformation
  "extract" "Extract meaningful information in prolate spheroidal coordinates"
  )
CreatePtkCommand(TensorImageToMesh
  "image2mesh" "Converts a tensor image into a tensor mesh structure"
  )
CreatePtkCommand(TensorMeshToImage
  "mesh2image" "Converts a tensor mesh structure into a tensor image"
  )
CreatePtkCommand(ExtrapolateTensorField
  "extrapolate" "Tensor Field Dense Approximation from Sparse Data"
  )
CreatePtkCommand(CreateSyntheticCardiacTensorMap
  "synthesize" "Create a synthetic tensor map out of some cardiac LV fibre structure information"
  )
CreatePtkCommand(FindClosestProlateTransform
  "find" "Find a prolate spheroid approximating an LV segmentation mesh"
  )
CreatePtkCommand(CreateProlateDomain
  "domain" "Create an ellipsoidal domain image out of a prolate sph. transform"
  )
CreatePtkCommand(TensorMeshStatistics
  "stats" "Compute statistics in prolate sph. coordinates over a set of tensors"
  )
CreatePtkCommand(CreateProlateAtlas
  "atlas" "Create an tensor atlas from a set of tensors in prolate coordinates"
  )
CreatePtkCommand(LimitImageToAHAZone
  "aha-image" "Crop an image to an AHA zone in the Prolate Spheroidal sense"
  )
CreatePtkCommand(LimitTensorImageToAHAZone
  "aha-tensors" "Crop a tensor image to an AHA zone in the Prolate Spheroidal sense"
  )

add_library(PTKCommands ${PTK_BUILD_TYPE}
  ${PTK_COMMAND_SRCS}
  )

target_link_libraries(PTKCommands
  ${ITK_LIBRARIES}
  PTKNumerics
  vtkCommon
  vtkIO
  ITKTensor
  )

if(NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
  install(TARGETS PTKCommands
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif(NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)