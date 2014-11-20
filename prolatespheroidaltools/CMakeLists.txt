project(vtkProlateTools)

if(WIN32)
  set(${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY lib)
  set(${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY bin)
  set(${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY bin)
else(WIN32)
  set(${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY lib)
  set(${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY bin)
  set(${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY lib)
endif(WIN32)

set(LIBRARY_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY})
set(ARCHIVE_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY})
set(RUNTIME_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})

set(LIBRARY_INSTALL_OUTPUT_PATH    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY})
set(ARCHIVE_INSTALL_OUTPUT_PATH    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY})
set(RUNTIME_INSTALL_OUTPUT_PATH    ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})
set(EXECUTABLE_INSTALL_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY})

find_package(ITK REQUIRED)
find_package(VTK REQUIRED)

if(ITK_FOUND AND VTK_FOUND)

  set(PROLATE_EXECUTABLES )
  
  include(${ITK_USE_FILE})
  include(${VTK_USE_FILE})

  set(${PROJECT_NAME}_INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
    )
  
  include_directories(
    ${${PROJECT_NAME}_INCLUDE_DIRS}
    )
  
  if (APPLE)
    set(GUI_EXECUTABLE MACOSX_BUNDLE)
  endif()

  add_library(vtkEllipsoidalSegmentation ${${PROJECT_NAME}_BUILD_TYPE}
    itkDoubleLandmarkSetToDomainImageFilter.cxx
    itkEllipsoidalTransformToDomainImageFilter.cxx
    itkOblateSpheroidalTransformToDomainImageFilter.cxx
    itkProlateSpheroidalTransformToDomainImageFilter.cxx
    vtkEllipsoidalTransformController.cxx
    vtkOblateSpheroidalTransformController.cxx
    vtkProlateSpheroidalTransformController.cxx
    )
  
  TARGET_LINK_LIBRARIES(vtkEllipsoidalSegmentation
    ${VTK_LIBRARIES}
    ${ITK_LIBRARIES}
    vtkLandmarkVariationalSegmentation
    )


  add_executable(prolategenerator ${GUI_EXECUTABLE}
    ProlateGenerator.cxx
    )
  target_link_libraries(prolategenerator
    vtkImageView
    vtkEllipsoidalSegmentation
    )
  set(PROLATE_EXECUTABLES prolategenerator ${PROLATE_EXECUTABLES} )
  
  add_executable(oblategenerator ${GUI_EXECUTABLE}
    OblateGenerator.cxx
    )
  target_link_libraries(oblategenerator
    vtkImageView
    vtkEllipsoidalSegmentation
    )
  set(PROLATE_EXECUTABLES oblategenerator ${PROLATE_EXECUTABLES} )
  add_executable(ellispoidalgenerator ${GUI_EXECUTABLE}
    EllipsoidalGenerator.cxx
    )
  target_link_libraries(ellispoidalgenerator
    vtkImageView
    vtkEllipsoidalSegmentation
    )
  set(PROLATE_EXECUTABLES ellispoidalgenerator ${PROLATE_EXECUTABLES} )

endif()

## #################################################################
## Installation
## #################################################################

if(NOT ${PROJECT_NAME}_INSTALL_NO_RUNTIME)
  foreach (exe ${PROLATE_EXECUTABLES})
    install(TARGETS ${exe}
      BUNDLE DESTINATION bin
      RUNTIME DESTINATION bin
      )
  endforeach()
  install_files(/include "(\\.h|\\.txx)$")
endif()

