#
# this module look for SPRNG (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# SPRNG_INCLUDE_PATH = where sprng.h can be found
# SPRNG_LIBRARY      = the library to link against (sprng etc)
#

SET(TRIAL_LIBRARY_PATHS
  /usr/local/lib
  /sw/lib
  $ENV{SPRNG_HOME}/lib
)
SET(TRIAL_INCLUDE_PATHS
  /usr/local/include
  /sw/include
  $ENV{SPRNG_HOME}/include
)

#IF($ENV{SPRNG_HOME} MATCHES "sprng")
#  SET(TRIAL_LIBRARY_PATHS
#     $ENV{SPRNG_HOME}/lib ${TRIAL_LIBRARY_PATHS} )
#  SET(TRIAL_INCLUDE_PATHS
#     $ENV{SPRNG_HOME}/include ${TRIAL_INCLUDE_PATHS} )
#ENDIF($ENV{SPRNG_HOME} MATCHES "sprng")

FIND_LIBRARY(SPRNG_LIBRARY
  sprng
  ${TRIAL_LIBRARY_PATHS}
)

FIND_PATH(SPRNG_INCLUDE_DIR
  sprng.h
  ${TRIAL_INCLUDE_PATHS} 
)

IF(SPRNG_INCLUDE_DIR AND SPRNG_LIBRARY)
  SET(FOUND_SPRNG 1 CACHE BOOL "Found sprng library")
ELSE(SPRNG_INCLUDE_DIR AND SPRNG_LIBRARY)
  SET(FOUND_SPRNG 0 CACHE BOOL "Not fount sprng library")
ENDIF(SPRNG_INCLUDE_DIR AND SPRNG_LIBRARY)

MARK_AS_ADVANCED(
  SPRNG_INCLUDE_DIR 
  SPRNG_LIBRARY 
  FOUND_SPRNG
)
