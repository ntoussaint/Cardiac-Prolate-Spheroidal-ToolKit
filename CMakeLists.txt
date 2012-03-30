project(CPSTK)


## #################################################################
## Required dependencies
## #################################################################

find_package(ITK REQUIRED)
find_package(VTK REQUIRED)
find_package(TTK REQUIRED)

include(${ITK_USE_FILE})
include(${VTK_USE_FILE})
include(${TTK_USE_FILE})


set(${PROJECT_NAME}_VERSION_MAJOR 1 CACHE STRING "CPSTK major version number.")
set(${PROJECT_NAME}_VERSION_MINOR 1 CACHE STRING "CPSTK minor version number.")
set(${PROJECT_NAME}_VERSION_BUILD 0 CACHE STRING "CPSTK build version number.")
set(${PROJECT_NAME}_VERSION
  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_BUILD}")

mark_as_advanced(${PROJECT_NAME}_VERSION_MAJOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_MINOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_BUILD)

## ###################################################################
## Some packages options are required for this project. Check them now
## ###################################################################

set(_ITKVersionString "${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH}" )
math(EXPR _ITKVersionNum "${ITK_VERSION_MAJOR}*100*100 + ${ITK_VERSION_MINOR}*100 + ${ITK_VERSION_PATCH}")
if( _ITKVersionNum LESS 30800 )
  message(SEND_ERROR "The ITK version you want to use (${_ITKVersionString}) is not supported by this project. Please use a more recent version of ITK. The minimum required version is 3.8.0. The recommended version is 3.20.1.")
endif( _ITKVersionNum LESS 30800 )
if(NOT ITK_USE_REVIEW)
  message(SEND_ERROR "ITK_USE_REVIEW has not been set to on when compiling ITK. This option is mandatory for this project to build. Please recompile your ITK installation.")
endif(NOT ITK_USE_REVIEW)

# -----------------------------------------------------------------------------
# Header file used to pass CMake settings to the source code
# -----------------------------------------------------------------------------
configure_file(${PROJECT_SOURCE_DIR}/cpstkConfigure.h.in
               ${PROJECT_BINARY_DIR}/cpstkConfigure.h CONFIGURE IMMEDIATE)

if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  install(FILES ${CPSTK_BINARY_DIR}/cpstkConfigure.h
      DESTINATION include)
endif(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)

#-----------------------------------------------------------------------------
# Add compiler flags CPSTK needs to work on this platform.  This must be
# done AFTER the call to CMAKE_EXPORT_BUILD_SETTINGS.
set(CPSTK_REQUIRED_C_FLAGS)
set(CPSTK_REQUIRED_CXX_FLAGS)
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${CPSTK_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CPSTK_REQUIRED_CXX_FLAGS}")

#-----------------------------------------------------------------------------
# Set any extra compilation flags here
if(CMAKE_COMPILER_IS_GNUCXX)
  # Add some useful warning flags
  # Note: -Wall and -Wno-deprecated are already pulled in by
  # UseITK.cmake
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wshadow")
endif(CMAKE_COMPILER_IS_GNUCXX)


include_directories(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/Numerics
  ${PROJECT_SOURCE_DIR}/Algorithms
  ${PROJECT_SOURCE_DIR}/Commands
  ${PROJECT_SOURCE_DIR}/IPF
)

include_directories(
  ${PROJECT_BINARY_DIR}
  ${PROJECT_BINARY_DIR}/Commands
)


## #################################################################
## BUILD TYPE
## #################################################################

set(CPSTK_BUILD_TYPE
STATIC
)

if (NOT WIN32)
  if (ITK_BUILD_SHARED)
    set (CPSTK_BUILD_TYPE SHARED)
  endif (ITK_BUILD_SHARED)
endif (NOT WIN32)


## #################################################################
## Actual CMakeLists code
## #################################################################

add_subdirectory(IPF)
add_subdirectory(Algorithms)
add_subdirectory(Numerics)
add_subdirectory(Commands)

add_executable(cpstk
  cpstk.cxx
  )

target_link_libraries(cpstk
  CPSTKCommands
  CPSTKProgramFactory
  ${ITK_LIBRARIES}
  )

## #################################################################
## Installation
## #################################################################

if(NOT ${PROJECT_NAME}_INSTALL_NO_RUNTIME)
  install(TARGETS cpstk
    BUNDLE DESTINATION bin
    RUNTIME DESTINATION bin
    )
endif(NOT ${PROJECT_NAME}_INSTALL_NO_RUNTIME)


#-----------------------------------------------------------------------------
# doxygen (and other?) documentation

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
