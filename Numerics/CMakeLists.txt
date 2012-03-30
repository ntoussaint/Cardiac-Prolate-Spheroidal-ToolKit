add_library(CPSTKNumerics ${CPSTK_BUILD_TYPE}
  itkKaiserBesselKernelFunction.cxx
  itkTensorMeshImageHybridCostFunction3.cxx
  itkNewUoaOptimizer.cxx
  )
target_link_libraries(CPSTKNumerics
  ${ITK_LIBRARIES}
  ITKTensor
  vtkIO
  )

if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  file(GLOB __files1 "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
  file(GLOB __files2 "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
  file(GLOB __files3 "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
  install(FILES ${__files1} ${__files2} ${__files3}
    DESTINATION include
    COMPONENT   Development
  )
endif(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
