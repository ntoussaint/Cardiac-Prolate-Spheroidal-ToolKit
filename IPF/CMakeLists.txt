add_library(ITKProgramFactory ${PTK_BUILD_TYPE}
itkCommandObjectBase.cxx
itkCommandObjectFactory.cxx
)

target_link_libraries(ITKProgramFactory
${ITK_LIBRARIES}
)

if (NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
  install(TARGETS ITKProgramFactory
  DESTINATION lib
  COMPONENT RuntimeLibraries
  )
endif (NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
