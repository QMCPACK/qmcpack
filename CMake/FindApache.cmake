#
# this module look for xerces-c from apache.org support
# it will define the following values
#
# XERCES_INCLUDE_DIR  = where xercesc/util/XMLString.hpp can be found
# XERCES_LIBRARY      = the library to link against xerces-c
# FOUND_XERCESC       = set to 1 if xerces-c is found
#
SET(TRIAL_LIBRARY_PATHS
  $ENV{XERCESC_HOME}/lib
  /usr/lib
  /usr/local/lib
  /sw/lib
  ) 
SET(TRIAL_INCLUDE_PATHS
  $ENV{XERCESC_HOME}/include
  /usr/include
  /usr/local/include
  /sw/include
  ) 

FIND_LIBRARY(XERCESC_LIBRARY xerces-c ${TRIAL_LIBRARY_PATHS})
FIND_PATH(XERCESC_INCLUDE_DIR xercesc/util/XMLString.hpp ${TRIAL_INCLUDE_PATHS})

IF(XERCESC_INCLUDE_DIR AND XERCESC_LIBRARY)
  SET(FOUND_XERCESC 1 CACHE BOOL "Found apache/xerces library")
ELSE(XERCESC_INCLUDE_DIR AND XERCESC_LIBRARY)
  SET(FOUND_XERCESC 0 CACHE BOOL "Missing apache/xerces library")
ENDIF(XERCESC_INCLUDE_DIR AND XERCESC_LIBRARY)

MARK_AS_ADVANCED(
  XERCESC_INCLUDE_DIR 
  XERCESC_LIBRARY 
  FOUND_XERCESC
  )

