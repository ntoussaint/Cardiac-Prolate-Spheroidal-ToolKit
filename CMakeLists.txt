cmake_minimum_required(VERSION 2.8)

set(CMAKE_PROJECT_NAME ${PROJECT_NAME})

project(CPSTK)

option(BUILD_GDCMIMPORT   "Build the GDCM DICOM Import tool"  OFF)
option(BUILD_MEDVIEW "Build a simple visualisation tool"      OFF)
option(BUILD_CPSTK "Build the CPSTK binaries"      OFF)
option(BUILD_SEGMENTATION "Build the Segmentation tool"       OFF)
option(BUILD_PROLATETOOLS "Build the Prolate Spheroidal Tools" OFF)

## #################################################################
## Configure path
## #################################################################


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

set(${PROJECT_NAME}_VERSION_MAJOR 1 CACHE STRING "${PROJECT_NAME} major version number.")
set(${PROJECT_NAME}_VERSION_MINOR 1 CACHE STRING "${PROJECT_NAME} minor version number.")
set(${PROJECT_NAME}_VERSION_BUILD 0 CACHE STRING "${PROJECT_NAME} build version number.")
set(${PROJECT_NAME}_VERSION
  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_BUILD}")

mark_as_advanced(${PROJECT_NAME}_VERSION_MAJOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_MINOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_BUILD)

# -----------------------------------------------------------------------------
# Header file used to pass CMake settings to the source code
# -----------------------------------------------------------------------------

configure_file(${PROJECT_SOURCE_DIR}/cpstkConfigure.h.in
               ${PROJECT_BINARY_DIR}/cpstkConfigure.h CONFIGURE IMMEDIATE)

if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  install(FILES ${${PROJECT_NAME}_BINARY_DIR}/cpstkConfigure.h
      DESTINATION include)
endif(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)

# To help using ${PROJECT_NAME} with other projects
configure_file(${${PROJECT_NAME}_SOURCE_DIR}/Use${PROJECT_NAME}.cmake.in
               ${${PROJECT_NAME}_BINARY_DIR}/Use${PROJECT_NAME}.cmake COPYONLY IMMEDIATE)

#-----------------------------------------------------------------------------
# Add compiler flags ${PROJECT_NAME} needs to work on this platform.  This must be
# done AFTER the call to CMAKE_EXPORT_BUILD_SETTINGS.
set(${PROJECT_NAME}_REQUIRED_C_FLAGS)
set(${PROJECT_NAME}_REQUIRED_CXX_FLAGS)
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${${PROJECT_NAME}_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${${PROJECT_NAME}_REQUIRED_CXX_FLAGS}")

#-----------------------------------------------------------------------------
# Set any extra compilation flags here
if(CMAKE_COMPILER_IS_GNUCXX)
  # Add some useful warning flags
  # Note: -Wall and -Wno-deprecated are already pulled in by
  # UseITK.cmake
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wshadow")
endif(CMAKE_COMPILER_IS_GNUCXX)



## #################################################################
## Required dependencies
## #################################################################

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

unset(${PROJECT_NAME}_THIRD_PARTY_PACKAGES)
set(${PROJECT_NAME}_THIRD_PARTY_PACKAGES
  VTK
  TTK
)

set(${PROJECT_NAME}_PACKAGES_MISSING 0)

foreach (package ${${PROJECT_NAME}_THIRD_PARTY_PACKAGES})
  find_package(${package} QUIET)
  if(${package}_FOUND)
    include(${${package}_USE_FILE})
  else(${package}_FOUND)
    set(${PROJECT_NAME}_PACKAGES_MISSING 1)
    set(${PROJECT_NAME}_MISSING_PACKAGES 
      ${package} 
      ${${PROJECT_NAME}_MISSING_PACKAGES}
      )
  endif(${package}_FOUND)
endforeach (package) 


#-----------------------------------------------------------------------------
# Set the include directories
# also include the binary directories as some source files are built there
#-----------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/Numerics
  ${PROJECT_SOURCE_DIR}/Algorithms
  ${PROJECT_SOURCE_DIR}/Commands
  ${PROJECT_SOURCE_DIR}/IPF
  ${PROJECT_SOURCE_DIR}/Utils
  ${PROJECT_SOURCE_DIR}/Utils/SymmetricElasticLogDemons
  ${PROJECT_BINARY_DIR}
  ${PROJECT_BINARY_DIR}/Commands
)

set(${CMAKE_PROJECT_NAME}_INCLUDE_DIRS
  ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
  ${PROJECT_NAME}_INCLUDE_DIRS
  )

include_directories(
  ${${PROJECT_NAME}_INCLUDE_DIRS}
)

#-----------------------------------------------------------------------------
# Include the external configuration file
include (${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake)


## #################################################################
## BUILD TYPE
## #################################################################

set(${PROJECT_NAME}_BUILD_TYPE
STATIC
)

if (NOT WIN32)
  if (ITK_BUILD_SHARED)
    set (${PROJECT_NAME}_BUILD_TYPE SHARED)
  endif (ITK_BUILD_SHARED)
endif (NOT WIN32)


## #################################################################
## Actual CMakeLists code
## #################################################################
  
if(NOT ${${PROJECT_NAME}_PACKAGES_MISSING})
  add_subdirectory(IPF)
  add_subdirectory(Algorithms)
  add_subdirectory(Numerics)
  add_subdirectory(Commands)
  
  add_executable(cpstk
    cpstk.cxx
    )
  target_link_libraries(cpstk
    ${PROJECT_NAME}Commands
    ${PROJECT_NAME}ProgramFactory
    ${ITK_LIBRARIES}
    )
else(NOT ${${PROJECT_NAME}_PACKAGES_MISSING})
  message(WARNING "Some packages are missing to completely build ${PROJECT_NAME} : ${${PROJECT_NAME}_MISSING_PACKAGES}. ${PROJECT_NAME} won't be built")  
endif(NOT ${${PROJECT_NAME}_PACKAGES_MISSING})



if (BUILD_MEDVIEW)
  add_subdirectory(vtkImageView)
  add_subdirectory(vtkDataManagement)

  set(${CMAKE_PROJECT_NAME}_INCLUDE_DIRS
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkImageView
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkDataManagement
    )

  add_subdirectory(medview)
endif()

if (BUILD_GDCMIMPORT)
  add_subdirectory(gdcmimport)
endif()

if (BUILD_SEGMENTATION)
  set(${CMAKE_PROJECT_NAME}_INCLUDE_DIRS
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/segmentation/LandmarkVariationalSegmentation
    )
  add_subdirectory(segmentation)
endif()

if (BUILD_PROLATETOOLS)
  add_subdirectory(prolatespheroidaltools)
endif()




## #################################################################
## Installation
## #################################################################

if(NOT ${PROJECT_NAME}_INSTALL_NO_RUNTIME)
  if (NOT ${${PROJECT_NAME}_PACKAGES_MISSING})
    install(TARGETS cpstk
      BUNDLE DESTINATION bin
      RUNTIME DESTINATION bin
      )
  endif()
endif()


#-----------------------------------------------------------------------------
# doxygen documentation

find_package(Doxygen)
if (DOXYGEN)
  option( LINK_EXTERNAL_DOC "Should the documentation be linked with external sources such as ITK?" NO )
  set(MY_LINK_EXTERNAL_DOC ${LINK_EXTERNAL_DOC})
  set(MY_DOXYGEN_BUILD_DIR ${PROJECT_BINARY_DIR}/Doxygen)
  
  
  configure_file (
    ${PROJECT_SOURCE_DIR}/Doxygen/doxygen.config.in
    ${MY_DOXYGEN_BUILD_DIR}/doxygen.config
  )

  if (MY_LINK_EXTERNAL_DOC)
    exec_program(${CMAKE_COMMAND} ${MY_DOXYGEN_BUILD_DIR}
      ARGS -E tar xvz ${PROJECT_SOURCE_DIR}/Doxygen/vtkNightlyDoc.tag.tar.gz
    )
    exec_program(${CMAKE_COMMAND} ${MY_DOXYGEN_BUILD_DIR}
      ARGS -E tar xvz ${PROJECT_SOURCE_DIR}/Doxygen/KWWidgetsNightlyDoc.tag.tar.gz
    )
    exec_program(${CMAKE_COMMAND} ${MY_DOXYGEN_BUILD_DIR}
      ARGS -E tar xvz ${PROJECT_SOURCE_DIR}/Doxygen/InsightDoxygen.tag.tar.gz
    )
  else (MY_LINK_EXTERNAL_DOC)
    exec_program(${CMAKE_COMMAND}
      ARGS -E remove ${MY_DOXYGEN_BUILD_DIR}/vtkNightlyDoc.tag
    )
    exec_program(${CMAKE_COMMAND}
      ARGS -E remove ${MY_DOXYGEN_BUILD_DIR}/KWWidgetsNightlyDoc.tag
    )
    exec_program(${CMAKE_COMMAND}
      ARGS -E remove ${MY_DOXYGEN_BUILD_DIR}/InsightDoxygen.tag
    )
  endif (MY_LINK_EXTERNAL_DOC)

  add_custom_target(doc
    ${DOXYGEN}
    ${MY_DOXYGEN_BUILD_DIR}/doxygen.config
  )
endif(DOXYGEN) 

add_subdirectory(Utils)
