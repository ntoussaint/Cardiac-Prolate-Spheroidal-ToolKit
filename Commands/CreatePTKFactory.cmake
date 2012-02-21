function(CreatePTKFactory src_dir command_name description)

  set(FILE_NAME  "itk${command_name}CommandFactory")
  set(DESCRIPTION ${description})
  string(TOUPPER ${FILE_NAME} HEADER_PROT)

  configure_file(${src_dir}/itkPtkCommandFactoryTemplate.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}.h @ONLY)
  configure_file(${src_dir}/itkPtkCommandFactoryTemplate.cxx.in
    ${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}.cxx @ONLY)

endfunction()

string(REPLACE " " ";" COMMAND_NAME ${MY_COMMAND_NAME})
string(REPLACE " " ";" DESCRIPTION  ${MY_DESCRIPTION})

CreatePTKFactory(${SRC_DIR} ${MY_COMMAND_NAME} ${MY_DESCRIPTION})
