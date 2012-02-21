project(PTK)


set(${PROJECT_NAME}_VERSION_MAJOR 1 CACHE STRING "PTK major version number.")
set(${PROJECT_NAME}_VERSION_MINOR 1 CACHE STRING "PTK minor version number.")
set(${PROJECT_NAME}_VERSION_BUILD 0 CACHE STRING "PTK build version number.")
set(${PROJECT_NAME}_VERSION
  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_BUILD}")

mark_as_advanced(${PROJECT_NAME}_VERSION_MAJOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_MINOR)
mark_as_advanced(${PROJECT_NAME}_VERSION_BUILD)


## #################################################################
## Required dependencies
## #################################################################

find_package(ITK REQUIRED)
find_package(VTK REQUIRED)
find_package(TTK REQUIRED)


## ###################################################################
## Some packages options are required for this project. Check them now
## ###################################################################

set(_ITKVersionString "${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH}" )
math(EXPR _ITKVersionNum "${ITK_VERSION_MAJOR}*100*100 + ${ITK_VERSION_MINOR}*100 + ${ITK_VERSION_PATCH}")
if( _ITKVersionNum LESS 30800 )
  message(SEND_ERROR "The ITK version you want to use (${_ITKVersionString}) is not supported by this project. Please use a more recent version of ITK. The minimum required version is 3.8.0. The recommended version is 3.14.0.")
endif( _ITKVersionNum LESS 30800 )
if(NOT ITK_USE_REVIEW)
  message(SEND_ERROR "ITK_USE_REVIEW has not been set to on when compiling ITK. This option is mandatory for this project to build. Please recompile your ITK installation.")
endif(NOT ITK_USE_REVIEW)



# -----------------------------------------------------------------------------
# Header file used to pass CMake settings to the source code
# -----------------------------------------------------------------------------
configure_file(${PROJECT_SOURCE_DIR}/ptkConfigure.h.in
               ${PROJECT_BINARY_DIR}/ptkConfigure.h CONFIGURE IMMEDIATE)

if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  install(FILES ${PTK_BINARY_DIR}/ptkConfigure.h
      DESTINATION include)
endif(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)

set(PTK_REQUIRED_C_FLAGS)
set(PTK_REQUIRED_CXX_FLAGS)

include(${ITK_USE_FILE})
include(${VTK_USE_FILE})
include(${TTK_USE_FILE})

include_directories(
  ${PROJECT_BINARY_DIR}
)
set(PTK_BUILD_TYPE
STATIC
)

if (NOT WIN32)
  if (ITK_BUILD_SHARED)
    set (PTK_BUILD_TYPE SHARED)
  endif (ITK_BUILD_SHARED)
endif (NOT WIN32)


#-----------------------------------------------------------------------------
# Set any extra compilation flags here
if(CMAKE_COMPILER_IS_GNUCXX)
  # Add some useful warning flags
  # Note: -Wall and -Wno-deprecated are already pulled in by
  # UseITK.cmake
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wshadow")
endif(CMAKE_COMPILER_IS_GNUCXX)

#-----------------------------------------------------------------------------
# Add compiler flags PTK needs to work on this platform.  This must be
# done AFTER the call to CMAKE_EXPORT_BUILD_SETTINGS.
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${PTK_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PTK_REQUIRED_CXX_FLAGS}")






add_subdirectory(IPF)
add_subdirectory(Algorithms)
add_subdirectory(Commands)



add_library(ITKProgramFactory ${PTK_BUILD_TYPE}
itkCommandObjectBase.cxx
itkCommandObjectFactory.cxx
)

target_link_libraries(ITKProgramFactory
${ITK_LIBRARIES}
)




set(PTK_COMMAND_SRCS
  itkTensorsToVTKCommandFactory.cxx
  itkRotateProlateSpheroidCommandFactory.cxx
  itkExtractProlateInformationCommandFactory.cxx
  itkTensorsToVTKCommand.cxx
  itkRotateProlateSpheroidCommand.cxx
  itkExtractProlateInformationCommand.cxx
  )

add_library(ptkCommands ${PTK_BUILD_TYPE}
  ${PTK_COMMAND_SRCS}
  )

target_link_libraries(ptkCommands
  ${ITK_LIBRARIES}
  vtkCommon
  vtkIO
  ITKTensor
  )





add_executable(ptk
  ptk.cxx
  )

target_link_libraries(ptk
  ptkCommands
  ITKProgramFactory
  ${ITK_LIBRARIES}
  )










if (NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)
  install(TARGETS ITKProgramFactory
  DESTINATION lib
  COMPONENT RuntimeLibraries
  )
endif (NOT ${PROJECT_NAME}_INSTALL_NO_LIBRARIES)


if(NOT ${PROJECT_NAME}_INSTALL_NO_RUNTIME)
  install(TARGETS ptk
    BUNDLE DESTINATION bin
    RUNTIME DESTINATION bin
    )
endif(NOT ${PROJECT_NAME}_INSTALL_NO_RUNTIME)



if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  file(GLOB __files1 "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
  file(GLOB __files2 "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
  file(GLOB __files3 "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
  install(FILES ${__files1} ${__files2} ${__files3}
    DESTINATION include
    COMPONENT   Development
  )
endif(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
