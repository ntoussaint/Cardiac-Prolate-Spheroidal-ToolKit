# Generate the CPSTKConfig.cmake file in the build tree.  Also configure
# one for installation.  The file tells external projects how to use
# CPSTK.

#-----------------------------------------------------------------------------
# Settings specific to the build tree.

# The install-only section is empty for the build tree.
SET(CPSTK_CONFIG_INSTALL_ONLY)

# The "use" file.
SET(CPSTK_USE_FILE ${CPSTK_BINARY_DIR}/UseCPSTK.cmake)

# The build settings file.
SET(CPSTK_BUILD_SETTINGS_FILE ${CPSTK_BINARY_DIR}/CPSTKBuildSettings.cmake)

# include directory.
set(CPSTK_INCLUDE_DIRS_CONFIG ${CPSTK_INCLUDE_DIRS})

# Library directory.
SET(CPSTK_LIBRARY_DIRS_CONFIG ${LIBRARY_OUTPUT_PATH})

# Runtime library directory.
SET(CPSTK_RUNTIME_LIBRARY_DIRS_CONFIG ${LIBRARY_OUTPUT_PATH})

#-----------------------------------------------------------------------------
# Configure CPSTKConfig.cmake for the build tree.
CONFIGURE_FILE(${CPSTK_SOURCE_DIR}/CPSTKConfig.cmake.in
               ${CPSTK_BINARY_DIR}/CPSTKConfig.cmake @ONLY IMMEDIATE)


#-----------------------------------------------------------------------------
# Settings specific to the install tree.

# The "use" file.
SET(CPSTK_USE_FILE ${CMAKE_INSTALL_PREFIX}/lib/UseCPSTK.cmake)

# The build settings file.
SET(CPSTK_BUILD_SETTINGS_FILE ${CMAKE_INSTALL_PREFIX}/lib/CPSTKBuildSettings.cmake)

# Link directories.
IF(CYGWIN AND BUILD_SHARED_LIBS)
  # In Cygwin programs directly link to the .dll files.
  SET(CPSTK_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/bin)
ELSE(CYGWIN AND BUILD_SHARED_LIBS)
  SET(CPSTK_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/lib)
ENDIF(CYGWIN AND BUILD_SHARED_LIBS)

# Runtime directories.
IF(WIN32)
  SET(CPSTK_RUNTIME_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/bin)
ELSE(WIN32)
  SET(CPSTK_RUNTIME_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}/lib)
ENDIF(WIN32)


#-----------------------------------------------------------------------------
# Configure CPSTKConfig.cmake for the install tree.
CONFIGURE_FILE(${CPSTK_SOURCE_DIR}/CPSTKConfig.cmake.in
               ${CPSTK_BINARY_DIR}/Utilities/CPSTKConfig.cmake @ONLY IMMEDIATE)


INSTALL_FILES(/lib FILES ${CPSTK_BINARY_DIR}/Utilities/CPSTKConfig.cmake)