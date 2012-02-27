add_library(PTKProgramFactory ${PTK_BUILD_TYPE}
itkCommandObjectBase.cxx
itkCommandObjectFactory.cxx
)

target_link_libraries(PTKProgramFactory
${ITK_LIBRARIES}
)

if (NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
  install(TARGETS PTKProgramFactory
  DESTINATION lib
  COMPONENT RuntimeLibraries
  )
endif (NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
