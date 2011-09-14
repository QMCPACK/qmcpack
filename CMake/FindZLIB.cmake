#
# Find the native ZLIB includes and library
#
# ZLIB_INCLUDE_DIR - where to find zlib.h, etc.
# ZLIB_LIBRARIES   - List of fully qualified libraries to link against when using zlib.
# ZLIB_FOUND       - Do not attempt to use zlib if "no" or undefined.

set(Libz z)
IF(QMC_BUILD_STATIC)
  set(Libz libz.a)
ENDIF(QMC_BUILD_STATIC)

#FIND_PATH(ZLIB_INCLUDE_DIR zlib.h ${ZLIB_HOME}/include $ENV{ZLIB_HOME}/include)
FIND_LIBRARY(ZLIB_LOC ${Libz} ${ZLIB_HOME}/lib64 ${ZLIB_HOME}/lib $ENV{ZLIB_HOME}/lib)

IF(ZLIB_LOC)
  set(ZLIB_LIBRARY ${ZLIB_LOC})
  set(ZLIB_FOUND TRUE)
ENDIF(ZLIB_LOC)

#IF(ZLIB_INCLUDE_DIR)
#  IF(ZLIB_LIBRARY)
#    SET( ZLIB_LIBRARIES ${ZLIB_LIBRARY} )
#    SET( ZLIB_FOUND TRUE)
#  ENDIF(ZLIB_LIBRARY)
#ENDIF(ZLIB_INCLUDE_DIR)
