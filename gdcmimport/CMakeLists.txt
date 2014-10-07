project(GDCMImport)

unset(${PROJECT_NAME}_THIRD_PARTY_PACKAGES)
set(${PROJECT_NAME}_THIRD_PARTY_PACKAGES
ITK
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

#-----------------------------------------------------------------------------
# Set the include directories
# also include the binary directories as some source files are built there
#-----------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}
)

include_directories(
  ${${PROJECT_NAME}_INCLUDE_DIRS}
)

## #################################################################
## BUILD TYPE
## #################################################################

set(${PROJECT_NAME}_BUILD_TYPE
SHARED
)
if (NOT WIN32)
  if (ITK_BUILD_SHARED)
    set (${PROJECT_NAME}_BUILD_TYPE SHARED)
  endif (ITK_BUILD_SHARED)
endif (NOT WIN32)

## #################################################################
## LIBRARIES
## #################################################################

set(SYSTEM_GDCM_LIBRARIES
  gdcmCommon
  gdcmDSED
  gdcmDICT
)

set(${PROJECT_NAME}_TARGET_LINK_LIBRARIES
  ${ITK_LIBRARIES}
  ${SYSTEM_GDCM_LIBRARIES}
)

if (WIN32 AND NOT CYGWIN AND NOT MINGW)
  add_library(getoptcompatlib STATIC
    getoptcompat/getopt.c
    getoptcompat/getopt_long.c
    )
  set(${PROJECT_NAME}_TARGET_LINK_LIBRARIES
    ${${PROJECT_NAME}_TARGET_LINK_LIBRARIES}
    getoptcompatlib
    )
endif()


if(${${PROJECT_NAME}_PACKAGES_MISSING})
  message(WARNING "  Some packages are missing to completely build ${PROJECT_NAME} : ${${PROJECT_NAME}_MISSING_PACKAGES}. ${PROJECT_NAME} won't be built")
else()
  if (NOT ${ITK_USE_SYSTEM_GDCM})
    message(FATAL_ERROR "  ITK needs to use SYSTEM GDCM to build gdcmimport")
  endif()
  
  ## #################################################################
  ## EXECUTABLE
  ## #################################################################
  
  add_executable(gdcmimport
    GDCMImport.cxx
    )
  target_link_libraries(gdcmimport
    ${${PROJECT_NAME}_TARGET_LINK_LIBRARIES}
    )
  
  ## #################################################################
  ## Installation
  ## #################################################################
  
  if(NOT ${PROJECT_NAME}_INSTALL_NO_RUNTIME)
    if (NOT ${${PROJECT_NAME}_PACKAGES_MISSING})
      install(TARGETS gdcmimport
	BUNDLE DESTINATION bin
	RUNTIME DESTINATION bin
	)
      install_files(/include "(\\.h|\\.txx)$")
    endif()
  endif()

endif()

