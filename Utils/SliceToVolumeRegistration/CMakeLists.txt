find_package(ITK REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)
add_executable(SliceToSliceRegistration
  itkSliceToSliceRegistration.cxx 
)
target_link_libraries(SliceToSliceRegistration
  ITKCommon 
  ITKIO
  ITKNumerics
  ITKStatistics 
)