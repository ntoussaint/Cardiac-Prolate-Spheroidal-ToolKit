add_library(CPSTKProgramFactory ${CPSTK_BUILD_TYPE}
itkCommandObjectBase.cxx
itkCommandObjectFactory.cxx
)

target_link_libraries(CPSTKProgramFactory
${ITK_LIBRARIES}
)

if (NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
  install(TARGETS CPSTKProgramFactory
  DESTINATION lib
  COMPONENT RuntimeLibraries
  )
endif (NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
