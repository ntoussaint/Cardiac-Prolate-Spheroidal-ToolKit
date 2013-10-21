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
)
