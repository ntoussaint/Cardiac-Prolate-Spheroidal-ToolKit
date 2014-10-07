
add_library(vtkLandmarkVariationalSegmentation ${${PROJECT_NAME}_BUILD_TYPE}
  vtkLandmarkWidget.cxx
  vtkLandmarkSegmentationController.cxx
  )

TARGET_LINK_LIBRARIES(vtkLandmarkVariationalSegmentation
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  )
