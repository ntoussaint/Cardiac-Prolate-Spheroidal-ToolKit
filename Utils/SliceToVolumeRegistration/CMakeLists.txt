find_package(ITK REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)
add_executable(slicetosliceregistration
  itkSliceToSliceRegistration.cxx 
)
target_link_libraries(slicetosliceregistration
  ${ITK_LIBRARIES}
)
add_executable(volumetovolumeregistration
  itkVolumeToVolumeRegistration.cxx 
)
target_link_libraries(volumetovolumeregistration
  ${ITK_LIBRARIES}
)
add_executable(slicetovolumeregistration
  itkSliceToVolumeRegistration.cxx
)
target_link_libraries(slicetovolumeregistration
  ${ITK_LIBRARIES}
  CPSTKNumerics
)


if(NOT ${PROJECT_NAME}_INSTALL_NO_RUNTIME)
  install(TARGETS slicetosliceregistration
    BUNDLE DESTINATION bin
    RUNTIME DESTINATION bin
    )
  install_files(/include "(\\.h|\\.txx)$")
  install(TARGETS volumetovolumeregistration
    BUNDLE DESTINATION bin
    RUNTIME DESTINATION bin
    )
  install_files(/include "(\\.h|\\.txx)$")
  install(TARGETS slicetovolumeregistration
    BUNDLE DESTINATION bin
    RUNTIME DESTINATION bin
    )
  install_files(/include "(\\.h|\\.txx)$")
endif()
