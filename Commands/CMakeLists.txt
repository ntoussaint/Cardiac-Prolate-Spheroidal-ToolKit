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
    ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.cxx
    PARENT_SCOPE
    )
  
endfunction()

set(PTK_COMMAND_SRCS
  itkTensorsToVTKCommand.cxx
  itkRotateProlateSpheroidCommand.cxx
  itkExtractProlateInformationCommand.cxx
  )

CreatePtkCommand(
  TensorsToVTK
  "itk2vtk" "Convert a tensor image or a list of tensor images into a vtkUnstructuredGrid structure"
  )
CreatePtkCommand(RotateProlateSpheroid
  "rotate" "Rotate a Prolate Spheroid according to a vtk file describing the antero-posterior line"
  )
CreatePtkCommand(ExtractProlateInformation
  "extract" "Extract meaningful information in prolate spheroidal coordinates"
  )

add_library(ptkCommands ${PTK_BUILD_TYPE}
  ${PTK_COMMAND_SRCS}
  )

target_link_libraries(ptkCommands
  ${ITK_LIBRARIES}
  vtkCommon
  vtkIO
  ITKTensor
  )

if(NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
  install(TARGETS ptkCommands
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif(NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)