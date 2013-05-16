find_package(ITK REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)
add_executable(slicetosliceregistration
  itkSliceToSliceRegistration.cxx 
)
target_link_libraries(slicetosliceregistration
  ITKCommon 
  ITKIO
  ITKNumerics
  ITKStatistics 
)
add_executable(volumetovolumeregistration
  itkVolumeToVolumeRegistration.cxx 
)
target_link_libraries(volumetovolumeregistration
  ITKCommon 
  ITKIO
  ITKNumerics
  ITKStatistics 
)
