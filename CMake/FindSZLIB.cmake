#
# Find the native SZLIB includes and library
#
# SZLIB_INCLUDE_DIR - where to find zlib.h, etc.
# SZLIB_LIBRARIES   - List of fully qualified libraries to link against when using zlib.
# SZLIB_FOUND       - Do not attempt to use zlib if "no" or undefined.

FIND_PATH(SZLIB_INCLUDE_DIR szlib.h
  $ENV{SZLIB_HOME}/include
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(SZLIB_LIBRARIES sz
  $ENV{SZLIB_HOME}/lib
  /usr/lib
  /usr/local/lib
)

SET(SZLIB_FOUND FALSE)
IF(SZLIB_LIBRARIES)
  SET(SZLIB_FOUND TRUE)
ENDIF(SZLIB_LIBRARIES)
