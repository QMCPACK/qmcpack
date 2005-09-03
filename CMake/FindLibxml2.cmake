#
# this module look for libxml (http://www.xmlsoft.org) support
# it will define the following values
#
# LIBXML2_INCLUDE_DIR  = where libxml/xpath.h can be found
# LIBXML2_LIBRARY      = the library to link against libxml2
# FOUND_LIBXML2        = set to 1 if libxml2 is found
#
IF(EXISTS ${PROJECT_CMAKE}/Libxml2Config.cmake)
  INCLUDE(${PROJECT_CMAKE}/Libxml2Config.cmake)
ENDIF(EXISTS ${PROJECT_CMAKE}/Libxml2Config.cmake)

IF(Libxml2_INCLUDE_DIRS)

  FIND_PATH(LIBXML2_INCLUDE_DIR libxml/xpath.h ${Libxml2_INCLUDE_DIRS})
  FIND_LIBRARY(LIBXML2_LIBRARY xml2 ${Libxml2_LIBRARY_DIRS})

ELSE(Libxml2_INCLUDE_DIRS)

  SET(TRIAL_LIBRARY_PATHS
    $ENV{LIBXML2_HOME}/lib
    /usr/lib
    /usr/local/lib
    /sw/lib
  ) 
  SET(TRIAL_INCLUDE_PATHS
    $ENV{LIBXML2_HOME}/include/libxml2
    /usr/include/libxml2
    /usr/local/include/libxml2
    /sw/include/libxml2
  ) 

  FIND_LIBRARY(LIBXML2_LIBRARY xml2 ${TRIAL_LIBRARY_PATHS})
  FIND_PATH(LIBXML2_INCLUDE_DIR libxml/xpath.h ${TRIAL_INCLUDE_PATHS})

ENDIF(Libxml2_INCLUDE_DIRS)

IF(LIBXML2_INCLUDE_DIR AND LIBXML2_LIBRARY)
  SET(FOUND_LIBXML2 1 CACHE BOOL "Found libxml2 library")
ELSE(LIBXML2_INCLUDE_DIR AND LIBXML2_LIBRARY)
  SET(FOUND_LIBXML2 0 CACHE BOOL "Not fount libxml2 library")
ENDIF(LIBXML2_INCLUDE_DIR AND LIBXML2_LIBRARY)

MARK_AS_ADVANCED(
  LIBXML2_INCLUDE_DIR 
  LIBXML2_LIBRARY 
  FOUND_LIBXML2
  )
