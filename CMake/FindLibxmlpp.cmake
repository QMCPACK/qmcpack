#
# this module look for libxml++ (http://libxmlplusplus.sourceforge.net) support
# it will define the following values
#
# LIBXMLPP_INCLUDE_PATH = where libxml++/libxml++.h can be found
# LIBXMLPP_LIBRARY      = the library to link against libxml++
# FOUND_LIBXMLPP        = set to 1 if libxml++ is found
#
SET(LIBXMLPP xml++-1.0)
IF(EXISTS ${PROJECT_CMAKE}/LibxmlppConfig.cmake)
  INCLUDE(${PROJECT_CMAKE}/LibxmlppConfig.cmake)
ENDIF(EXISTS ${PROJECT_CMAKE}/LibxmlppConfig.cmake)

IF(Libxmlpp_INCLUDE_DIRS)

  FIND_PATH(LIBXMLPP_INCLUDE_DIR libxml++/libxml++.h ${Libxmlpp_INCLUDE_DIRS})
  FIND_LIBRARY(LIBXMLPP_LIBRARY ${LIBXMLPP} ${Libxmlpp_LIBRARY_DIRS})

ELSE(Libxmlpp_INCLUDE_DIRS)

  FIND_LIBRARY(LIBXMLPP_LIBRARY 
    ${LIBXMLPP}
    /usr/app/lib 
    /usr/lib 
    /usr/local/lib
    /sw/lib
    )

  FIND_PATH(LIBXMLPP_INCLUDE_DIR 
    libxml++/libxml++.h
    /usr/app/include/libxml++-1.0
    /usr/include/libxml++-1.0
    /usr/local/include/libxml++-1.0
    /sw/include/libxml++-1.0
    )
ENDIF(Libxmlpp_INCLUDE_DIRS)

IF(LIBXMLPP_INCLUDE_DIR AND LIBXMLPP_LIBRARY)
  SET(FOUND_LIBXMLPP 1 CACHE BOOL "Found libxml++ library")
ELSE(LIBXMLPP_INCLUDE_DIR AND LIBXMLPP_LIBRARY)
  SET(FOUND_LIBXMLPP 0 CACHE BOOL "Not fount libxml2 library")
ENDIF(LIBXMLPP_INCLUDE_DIR AND LIBXMLPP_LIBRARY)

MARK_AS_ADVANCED(
  LIBXMLPP_INCLUDE_DIR 
  LIBXMLPP_LIBRARY 
  FOUND_LIBXMLPP
  )
