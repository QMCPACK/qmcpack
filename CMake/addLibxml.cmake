#
# this module look for libxml (http://http://xmlsoft.org) support
# it will define the following values
#
# LIBXML_INCLUDE_PATH = where libxml/parser.h can be found
# LIBXML_LIBRARY      = the library to link against libxml2
#

FIND_LIBRARY(LIBXML_LIBRARY xml2 
             /usr/lib 
             /usr/local/lib
             /sw/lib)
FIND_PATH(LIBXML_INCLUDE_DIR 
          libxml/parser.h 
          /usr/include/libxml2
          /sw/include/libxml2
)

IF(NOT LIBXML_LIBRARY)
  MESSAGE(STATUS "libxml is not found in pre-defined paths.")
  IF(UNIX)
    IF(RUN_CONFIGURE)
      EXEC_PROGRAM(
       "${AUXPACAKGES}/libxml2/configure --prefix=${PROJECT_BINARY_DIR} --with-cc=${CMAKE _CC_COMPILER}"
        ${AUXPACAKGES}/libxml2
      )
      EXEC_PROGRAM(make ${AUXPACKAGES}/libxml2)
      EXEC_PROGRAM(make ${AUXPACKAGES}/libxml2 ARGS install)
      FIND_LIBRARY(LIBXML_LIBRARY xml2 ${PROJECT_BINARY_DIR}/lib)
      FIND_PATH(LIBXML_INCLUDE_DIR libxml/parser.h ${PROJECT_BINARY_DIR}/include)
    ENDIF(RUN_CONFIGURE)
  ENDIF(UNIX)
ENDIF(NOT LIBXML_LIBRARY)

IF(LIBXML_LIBRARY)
  MARK_AS_ADVANCED(LIBXML_INCLUDE_DIR LIBXML_LIBRARY)
  ADD_DEFINITIONS(-DENABLE_LIBXML2)
  INCLUDE_DIRECTORIES(${LIBXML_INCLUDE_DIR})
  LINK_LIBRARIES(${LIBXML_LIBRARY})
ENDIF(LIBXML_LIBRARY)
