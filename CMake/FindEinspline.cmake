#
# this module look for einspline support
# it will define the following values
#
# EINSPLINE_INCLUDE_DIR  = where einspline/config.h can be found
# EINSPLINE_LIBRARIES    = the library to link against 
# EINSPLINE_FOUND        = set to true after finding the pre-compiled library
#
# EINSPLINE_DIR          = directory where the source files are located
#
#FIND_LIBRARY(EINSPLINE_LIBRARY einspline ${TRIAL_LIBRARY_PATHS})
#FIND_PATH(EINSPLINE_INCLUDE_DIR einspline/bspline.h ${TRIAL_INCLUDE_PATHS} )
FIND_LIBRARY(EINSPLINE_LIBRARIES einspline ${QMC_LIBRARY_PATHS})
FIND_PATH(EINSPLINE_INCLUDE_DIR einspline/bspline.h ${QMC_INCLUDE_PATHS})

SET(EINSPLINE_FOUND FALSE)

IF(EINSPLINE_INCLUDE_DIR AND EINSPLINE_LIBRARIES)
  SET(EINSPLINE_FOUND TRUE)
  SET(HAVE_EINSPLINE_EXT 1)
  MARK_AS_ADVANCED(
    EINSPLINE_INCLUDE_DIR 
    EINSPLINE_LIBRARIES 
    EINSPLINE_FOUND
    )
ELSE()
  MESSAGE("-- Cannot find einspline library. Try to compile with the sources")
  FIND_PATH(EINSPLINE_SRC_DIR src/bspline.h ${EINSPLINE_HOME} $ENV{EINSPLINE_HOME})
  IF(EINSPLINE_SRC_DIR)
    SET(EINSPLINE_FOUND TRUE)
    SET(EINSPLINE_DIR ${EINSPLINE_SRC_DIR})
    SET(EINSPLINE_INCLUDE_DIR ${EINSPLINE_SRC_DIR})
    MARK_AS_ADVANCED(
      EINSPLINE_INCLUDE_DIR 
      EINSPLINE_DIR 
      EINSPLINE_FOUND
     )
  ELSE(EINSPLINE_SRC_DIR)
    MESSAGE("-- EINSPLINE_HOME is not found. Disable Einspline library.")
    MESSAGE("-- Download einspline library to utilize an optimized 3D-bspline library.")
  ENDIF(EINSPLINE_SRC_DIR)
ENDIF()

