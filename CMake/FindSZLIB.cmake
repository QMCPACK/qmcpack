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

FIND_LIBRARY(SZLIB_LIBRARY sz
  $ENV{SZLIB_HOME}/lib
  /usr/lib
  /usr/local/lib
)

IF(SZLIB_INCLUDE_DIR)
  IF(SZLIB_LIBRARY)
    SET( SZLIB_FOUND "YES" )
  ENDIF(SZLIB_LIBRARY)
ENDIF(SZLIB_INCLUDE_DIR)
