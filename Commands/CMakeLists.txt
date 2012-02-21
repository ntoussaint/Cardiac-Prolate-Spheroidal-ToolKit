function(CreateCommandFactory command_name description src_list)
  set(FACTORY_NAME  "itk${command_name}CommandFactory")
  
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.h ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.cxx
    COMMAND ${CMAKE_COMMAND} 
    -D SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR} 
    -D MY_COMMAND_NAME=${command_name}
    -D MY_DESCRIPTION=${description} 
    -P ${CMAKE_CURRENT_SOURCE_DIR}/CreatePTKFactory.cmake
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/itkPtkCommandFactoryTemplate.h.in
            ${CMAKE_CURRENT_SOURCE_DIR}/itkPtkCommandFactoryTemplate.h.in)
	  
  set(FACTORY_SRC_NAMES "${FACTORY_SRC_NAMES}; ${FACTORY_NAME}")
  set(${src_list}
    ${${src_list}}
    ${CMAKE_CURRENT_BINARY_DIR}/${FACTORY_NAME}.cxx
    PARENT_SCOPE
    )
  message("ptk_command_srcs:"${${src_list}})
  
endfunction()

set(FACTORY_SRC_NAMES )

set(PTK_COMMAND_SRCS
  itkTensorsToVTKCommand.cxx
  itkRotateProlateSpheroidCommand.cxx
  itkExtractProlateInformationCommand.cxx
  )

CreateCommandFactory(TensorsToVTK    "analyze"           PTK_COMMAND_SRCS)
CreateCommandFactory(RotateProlateSpheroid "analyze"     PTK_COMMAND_SRCS)
CreateCommandFactory(ExtractProlateInformation "analyze" PTK_COMMAND_SRCS)



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