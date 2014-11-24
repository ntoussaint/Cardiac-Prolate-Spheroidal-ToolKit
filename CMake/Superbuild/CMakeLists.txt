project(CPSTK_SUPERBUILD)
cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

option(BUILD_SHARED_LIBS "To build the libraries as shared" OFF)

set(CMAKE_PROPAGATED_VARIABLES
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
  -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  )

#-----------------------------------------------------------------------------
# Search for VTK. Download and install if not found
#-----------------------------------------------------------------------------
option(USE_SYSTEM_VTK "Use an already installed version of VTK" OFF)
if(USE_SYSTEM_VTK)
  find_package(VTK REQUIRED)
else()
  ExternalProject_Add(vtk
    GIT_REPOSITORY http://vtk.org/VTK.git
    GIT_TAG "v6.1.0"
    PREFIX ${PROJECT_BINARY_DIR}/vtk
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/vtk/download
    SOURCE_DIR ${CMAKE_BINARY_DIR}/vtk/source
    STAMP_DIR ${PROJECT_BINARY_DIR}/vtk/stamps
    TMP_DIR ${PROJECT_BINARY_DIR}/vtk/tmp
    BINARY_DIR ${PROJECT_BINARY_DIR}/vtk/build
    TIMEOUT 100
    UPDATE_COMMAND ""
    CMAKE_ARGS 
      ${CMAKE_PROPAGATED_VARIABLES} 
      -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/vtk/install
      -DBUILD_TESTING=OFF
      -DBUILD_EXAMPLES=OFF
      -DVTK_REQUIRED_OBJCXX_FLAGS=
    )
  set(VTK_DIR ${PROJECT_BINARY_DIR}/vtk/build)
  set(VTK_INSTALL_DIR ${PROJECT_BINARY_DIR}/vtk/install)
  install(DIRECTORY ${VTK_INSTALL_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(DIRECTORY ${VTK_INSTALL_DIR}/lib DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(DIRECTORY ${VTK_INSTALL_DIR}/bin DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
endif()



#-----------------------------------------------------------------------------
# Search for GDCM. Download and install if not found
#-----------------------------------------------------------------------------
option(USE_SYSTEM_GDCM "Use an already installed version of GDCM" OFF)
if(USE_SYSTEM_GDCM)
  find_package(GDCM REQUIRED)
else()
  ExternalProject_Add(gdcm
    GIT_REPOSITORY http://git.code.sf.net/p/gdcm/gdcm
    GIT_TAG "v2.4.4"
    PREFIX ${PROJECT_BINARY_DIR}/gdcm
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/gdcm/download
    SOURCE_DIR ${CMAKE_BINARY_DIR}/gdcm/source
    STAMP_DIR ${PROJECT_BINARY_DIR}/gdcm/stamps
    TMP_DIR ${PROJECT_BINARY_DIR}/gdcm/tmp
    BINARY_DIR ${PROJECT_BINARY_DIR}/gdcm/build
    TIMEOUT 100
    UPDATE_COMMAND ""
    CMAKE_ARGS 
      ${CMAKE_PROPAGATED_VARIABLES} 
      -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/gdcm/install
      -DBUILD_TESTING=OFF
      -DGDCM_BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
    )
  set(GDCM_DIR ${PROJECT_BINARY_DIR}/gdcm/build)
  set(GDCM_INSTALL_DIR ${PROJECT_BINARY_DIR}/gdcm/install)
  install(DIRECTORY ${GDCM_INSTALL_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(DIRECTORY ${GDCM_INSTALL_DIR}/lib DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()


#-----------------------------------------------------------------------------
# Search for ITK. Download and install if not found
#-----------------------------------------------------------------------------
option(USE_SYSTEM_ITK "Use an already installed version of ITK" OFF)
if(USE_SYSTEM_ITK)
  find_package(ITK REQUIRED)
else()
  ExternalProject_Add(itk
    GIT_REPOSITORY http://itk.org/ITK.git
    GIT_TAG "v4.6.1"
    PREFIX ${PROJECT_BINARY_DIR}/itk
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/itk/download
    SOURCE_DIR ${CMAKE_BINARY_DIR}/itk/source
    STAMP_DIR ${PROJECT_BINARY_DIR}/itk/stamps
    TMP_DIR ${PROJECT_BINARY_DIR}/itk/tmp
    BINARY_DIR ${PROJECT_BINARY_DIR}/itk/build
    TIMEOUT 100
    UPDATE_COMMAND ""
    CMAKE_ARGS 
      ${CMAKE_PROPAGATED_VARIABLES} 
      -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/itk/install
      -DITK_USE_SYSTEM_GDCM=ON
      -DGDCM_DIR=${GDCM_DIR}
      -DModule_ITKVtkGlue=ON
      -DVTK_DIR=${VTK_DIR}
      -DModule_ITKIOPhilipsREC=ON
      -DModule_ITKReview=ON
      -DBUILD_TESTING=OFF
      -DBUILD_EXAMPLES=OFF
    )
  add_dependencies(itk vtk gdcm)
  set(ITK_DIR ${PROJECT_BINARY_DIR}/itk/build)
  set(ITK_INSTALL_DIR ${PROJECT_BINARY_DIR}/itk/install)
  install(DIRECTORY ${ITK_INSTALL_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(DIRECTORY ${ITK_INSTALL_DIR}/lib DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(DIRECTORY ${ITK_INSTALL_DIR}/bin DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
endif()



#-----------------------------------------------------------------------------
# Search for TTK. Download and install if not found
#-----------------------------------------------------------------------------
option(USE_SYSTEM_TTK "Use an already installed version of TTK" OFF)
if(USE_SYSTEM_TTK)
  find_package(TTK REQUIRED)
else()
  ExternalProject_Add(ttk
    GIT_REPOSITORY https://github.com/ntoussaint/TTK-Public.git
    GIT_TAG "origin/experimental-ntoussaint"
    PREFIX ${PROJECT_BINARY_DIR}/ttk
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/ttk/download
    SOURCE_DIR ${CMAKE_BINARY_DIR}/ttk/source
    STAMP_DIR ${PROJECT_BINARY_DIR}/ttk/stamps
    TMP_DIR ${PROJECT_BINARY_DIR}/ttk/tmp
    BINARY_DIR ${PROJECT_BINARY_DIR}/ttk/build
    TIMEOUT 100
    UPDATE_COMMAND ""
    CMAKE_ARGS 
      ${CMAKE_PROPAGATED_VARIABLES} 
      -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/ttk/install
      -DBUILD_TESTING=OFF
      -DITK_DIR=${ITK_DIR}
      -DVTK_DIR=${VTK_DIR}
    )
  add_dependencies(ttk itk vtk)
  set(TTK_DIR ${PROJECT_BINARY_DIR}/ttk/build)
  set(TTK_INSTALL_DIR ${PROJECT_BINARY_DIR}/ttk/install)
  install(DIRECTORY ${TTK_INSTALL_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(DIRECTORY ${TTK_INSTALL_DIR}/lib DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(DIRECTORY ${TTK_INSTALL_DIR}/bin DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
endif()



#-----------------------------------------------------------------------------
# Search for CPSTK. Download and install if not found
#-----------------------------------------------------------------------------
option(USE_SYSTEM_CPSTK "Use an already installed version of CPSTK" OFF)
if(USE_SYSTEM_CPSTK)
  find_package(CPSTK REQUIRED)
else()
  ExternalProject_Add(cpstk
    GIT_REPOSITORY https://github.com/ntoussaint/Cardiac-Prolate-Spheroidal-ToolKit.git
    GIT_TAG "origin/master"
    PREFIX ${PROJECT_BINARY_DIR}/cpstk
    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/cpstk/download
    SOURCE_DIR ${CMAKE_BINARY_DIR}/cpstk/source
    STAMP_DIR ${PROJECT_BINARY_DIR}/cpstk/stamps
    TMP_DIR ${PROJECT_BINARY_DIR}/cpstk/tmp
    BINARY_DIR ${PROJECT_BINARY_DIR}/cpstk/build
    TIMEOUT 100
    UPDATE_COMMAND ""
    CMAKE_ARGS 
      ${CMAKE_PROPAGATED_VARIABLES} 
      -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/cpstk/install
      -DBUILD_TESTING=OFF
      -DITK_DIR=${ITK_DIR}
      -DVTK_DIR=${VTK_DIR}
      -DTTK_DIR=${TTK_DIR}
      -DBUILD_CPSTK=ON
      -DBUILD_GDCMIMPORT=ON
      -DBUILD_MEDVIEW=ON
      -DBUILD_PROLATETOOLS=ON
      -DBUILD_SEGMENTATION=ON
    )
  add_dependencies(cpstk itk vtk ttk)
  set(CPSTK_DIR ${PROJECT_BINARY_DIR}/cpstk/build)
  set(CPSTK_INSTALL_DIR ${PROJECT_BINARY_DIR}/cpstk/install)
  install(DIRECTORY ${CPSTK_INSTALL_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(DIRECTORY ${CPSTK_INSTALL_DIR}/lib DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(DIRECTORY ${CPSTK_INSTALL_DIR}/bin DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
endif()


