function(CreateCommandFactory command_name keyword description src_list)
  set(COMMAND_NAME  "itk${command_name}Command")
  set(FACTORY_NAME  "itk${command_name}CommandFactory")
  
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${COMMAND_NAME}.h ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.h ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.cxx
    COMMAND ${CMAKE_COMMAND} 
    -D SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR} 
    -D MY_COMMAND_NAME=${command_name}
    -D MY_DESCRIPTION=${description} 
    -D MY_KEYWORD=${keyword}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/CreatePTKCommand.cmake
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/itkPtkCommand.h.in
	    ${CMAKE_CURRENT_SOURCE_DIR}/itkPtkCommandFactoryTemplate.h.in
            ${CMAKE_CURRENT_SOURCE_DIR}/itkPtkCommandFactoryTemplate.cxx.in
	    )
	  
  set(${src_list}
    ${${src_list}}
    ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.cxx
    PARENT_SCOPE
    )
  
endfunction()

set(PTK_COMMAND_SRCS
  itkTensorsToVTKCommand.cxx
  itkRotateProlateSpheroidCommand.cxx
  itkExtractProlateInformationCommand.cxx
  )

CreatePTKCommand(TensorsToVTK
  "itk2vtk" ""
  PTK_COMMAND_SRCS
  )
CreatePTKCommand(RotateProlateSpheroid
  "rotate" ""
  PTK_COMMAND_SRCS
  )
CreatePTKCommand(ExtractProlateInformation
  "extract" ""
  PTK_COMMAND_SRCS
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