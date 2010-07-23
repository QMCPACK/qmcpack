#
# this module look for libxml (http://www.xmlsoft.org) support
# it will define the following values
#
# LibXml2_INCLUDE_DIR  = where libxml/xpath.h can be found
# LibXml2_LIBRARIES    = the library to link against libxml2
# LibXml2_FOUND        = set to 1 if libxml2 is found
#
IF(EXISTS ${PROJECT_CMAKE}/Libxml2Config.cmake)
  INCLUDE(${PROJECT_CMAKE}/Libxml2Config.cmake)
ENDIF(EXISTS ${PROJECT_CMAKE}/Libxml2Config.cmake)

IF(Libxml2_INCLUDE_DIRS)

  FIND_PATH(LIBXML2_INCLUDE_DIR libxml/xpath.h ${Libxml2_INCLUDE_DIRS})
  FIND_LIBRARY(LIBXML2_LIBRARY xml2 ${Libxml2_LIBRARY_DIRS})

ELSE(Libxml2_INCLUDE_DIRS)

  FIND_PATH(LIBXML_INCLUDE_DIR libxml2/libxml/xpath.h PATHS ${QMC_INCLUDE_PATHS} /usr/local/include /usr/include NO_DEFAULT_PATH)
  if(LIBXML_INCLUDE_DIR)
    set(LIBXML2_INCLUDE_DIR "${LIBXML_INCLUDE_DIR}/libxml2")
  else(LIBXML_INCLUDE_DIR)
    FIND_PATH(LIBXML2_INCLUDE_DIR libxml/xpath.h ${QMC_INCLUDE_PATHS})
  endif(LIBXML_INCLUDE_DIR)

  if(CPU_FLAGS MATCHES "lm")
    FIND_LIBRARY(LIBXML2_LIBRARIES xml2 PATHS ${QMC_LIBRARY_PATHS} /usr/local/lib64 /usr/lib64 /usr/local/lib /usr/lib NO_DEFAULT_PATH)
  else(CPU_FLAGS MATCHES "lm")
    FIND_LIBRARY(LIBXML2_LIBRARIES xml2 PATHS ${QMC_LIBRARY_PATHS} /usr/local/lib /usr/lib NO_DEFAULT_PATH)
  endif(CPU_FLAGS MATCHES "lm")

ENDIF(Libxml2_INCLUDE_DIRS)

SET(LIBXML2_FOUND FALSE)
IF(LIBXML2_INCLUDE_DIR AND LIBXML2_LIBRARIES)
  message(STATUS "LIBXML2_INCLUDE_DIR="${LIBXML2_INCLUDE_DIR})
  message(STATUS "LIBXML2_LIBRARIES="${LIBXML2_LIBRARIES})
  SET(LIBXML2_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
  LIBXML2_INCLUDE_DIR 
  LIBXML2_LIBRARIES 
  LIBXML2_FOUND
  )
