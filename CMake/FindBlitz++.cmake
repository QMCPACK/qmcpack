#
# this module look for blitz++ (http://www.oonumerics.org/blitz) support
# it will define the following values
#
# BLITZ_INCLUDE_DIR = where blitz/blitz.h can be found
#
# May want to define this but seldom required
# BLITZ_LIBRARY = where blitz library can be found (reserved)
#
#IF(EXISTS ${PROJECT_CMAKE}/BlitzppConfig.cmake)
#  INCLUDE(${PROJECT_CMAKE}/BlitzppConfig.cmake)
#ENDIF(EXISTS ${PROJECT_CMAKE}/BlitzppConfig.cmake)

IF(Blitzpp_INCLUDE_DIRS)
  FIND_PATH(BLITZ_INCLUDE_DIR blitz/blitz.h  ${Blitzpp_INCLUDE_DIRS})
ELSE(Blitzpp_INCLUDE_DIRS)
  SET(TRIAL_PATHS
    $ENV{BLITZ_HOME}
    /usr/apps/include
    /usr/include
    /opt/include
    /usr/local/include
  )

  FIND_PATH(BLITZ_INCLUDE_DIR blitz/blitz.h ${TRIAL_PATHS})
ENDIF(Blitzpp_INCLUDE_DIRS)

IF(BLITZ_INCLUDE_DIR)
  SET(FOUND_BLITZ 1 CACHE BOOL "Found blitz++ library")
ELSE(BLITZ_INCLUDE_DIR)
  SET(FOUND_BLITZ 0 CACHE BOOL "Found blitz++ library")
ENDIF(BLITZ_INCLUDE_DIR)

MARK_AS_ADVANCED(
   BLITZ_INCLUDE_DIR
   FOUND_BLITZ
)
