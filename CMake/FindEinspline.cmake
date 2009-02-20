#
# this module look for einspline support
# it will define the following values
#
# EINSPLINE_INCLUDE_DIR  = where einspline/config.h can be found
# EINSPLINE_LIBRARY      = the library to link against 
# FOUND_EINSPLINE        = set to true after finding the library
#
SET(TRIAL_LIBRARY_PATHS
    $ENV{EINSPLINE_HOME}/lib
    /usr/apps/lib
    /usr/lib 
    /usr/local/lib
    /opt/lib
    /sw/lib
   )

SET(TRIAL_INCLUDE_PATHS
    $ENV{EINSPLINE_HOME}/include
    /usr/apps/include
    /usr/include
    /opt/include
    /usr/local/include
    /sw/include
   )

FIND_LIBRARY(EINSPLINE_LIBRARY einspline ${TRIAL_LIBRARY_PATHS})
FIND_PATH(EINSPLINE_INCLUDE_DIR einspline/bspline.h ${TRIAL_INCLUDE_PATHS} )

IF(EINSPLINE_INCLUDE_DIR AND EINSPLINE_LIBRARY)
  SET(FOUND_EINSPLINE 1 CACHE BOOL "Found einspline pre-buidlt library")
  SET(HAVE_EINSPLINE 1)
  SET(HAVE_EINSPLINE_EXT 1)
  MARK_AS_ADVANCED(
      EINSPLINE_INCLUDE_DIR 
      EINSPLINE_LIBRARY 
      FOUND_EINSPLINE
      )
ELSE(EINSPLINE_INCLUDE_DIR AND EINSPLINE_LIBRARY)
  FIND_PATH(EINSPLINE_SRC_DIR src/bspline.h ${TRIAL_INCLUDE_PATHS} )
  IF(EINSPLINE_SRC_DIR)
  ELSE(EINSPLINE_SRC_DIR)
    SET(FOUND_EINSPLINE 0 CACHE BOOL "Cannot find einspline source or library.")
    MESSAGE("-- EINSPLINE_HOME is not found. Disable Einspline library.")
    MESSAGE("-- Download einspline library to utilize an optimized 3D-bspline library.")
    SET(HAVE_EINSPLINE 1)
  ENDIF(EINSPLINE_SRC_DIR)
ENDIF(EINSPLINE_INCLUDE_DIR AND EINSPLINE_LIBRARY)

