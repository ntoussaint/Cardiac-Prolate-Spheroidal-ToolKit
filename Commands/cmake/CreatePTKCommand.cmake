function(CreatePTKCommand command_name keyword description src_dir)

  set(FILE_NAME  "itk${command_name}CommandFactory")
  set(DESCRIPTION ${description})
  string(TOUPPER ${FILE_NAME} HEADER_PROT)

  configure_file(${src_dir}/itkPtkCommandFactoryTemplate.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}.h @ONLY)
  configure_file(${src_dir}/itkPtkCommandFactoryTemplate.cxx.in
    ${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}.cxx @ONLY)
  
  set(FILE_NAME  "itk${command_name}Command")
  string(TOUPPER ${FILE_NAME} HEADER_PROT)
  configure_file(${src_dir}/itkPtkCommandTemplate.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME}.h @ONLY)
  
endfunction()

string(REPLACE " " ";" COMMAND_NAME ${MY_COMMAND_NAME})
string(REPLACE " " ";" DESCRIPTION  ${MY_DESCRIPTION})
string(REPLACE " " ";" KEYWORD  ${MY_KEYWORD})

CreatePTKCommand(${MY_COMMAND_NAME} ${MY_KEYWORD} ${MY_DESCRIPTION} ${SRC_DIR})
