project(vtkSegmentationToolKit)

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

  include(${ITK_USE_FILE})
  include(${VTK_USE_FILE})

  set(${PROJECT_NAME}_INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/LandmarkVariationalSegmentation
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
    )

  include_directories(
    ${${PROJECT_NAME}_INCLUDE_DIRS}
    )
  
  add_subdirectory(LandmarkVariationalSegmentation)
  
  if (APPLE)
    set(GUI_EXECUTABLE MACOSX_BUNDLE)
  endif()

  add_executable(segmentationtool ${GUI_EXECUTABLE}
    SegmentationTool.cxx
    )
  target_link_libraries(segmentationtool
    vtkImageView
    vtkLandmarkVariationalSegmentation
    )


endif()

## #################################################################
## Installation
## #################################################################

if(NOT ${PROJECT_NAME}_INSTALL_NO_RUNTIME)
  install(TARGETS segmentationtool
    BUNDLE DESTINATION bin
    RUNTIME DESTINATION bin
    )
  install_files(/include "(\\.h|\\.txx)$")
endif()

