#
# this module look for libxml++ (http://libxmlplusplus.sourceforge.net) support
# it will define the following values
#
# LIBXMLPP_INCLUDE_PATH = where libxml++/libxml++.h can be found
# LIBXMLPP_LIBRARY      = the library to link against
#

FIND_LIBRARY(LIBXMLPP_LIBRARY xml++-1.0 
             lib
             /usr/lib 
             /usr/local/lib
             /usr/apps/lib
             /usr/apps/tools/libxml++/lib
             )
FIND_PATH(LIBXMLPP_INCLUDE_DIR
         libxml++/libxml++.h
         /usr/apps/include/libxml++-1.0
         /usr/local/include/libxml++-1.0
         /usr/apps/tools/libxml++/include
         libxml++/libxml++.h
         include/libxml++-1.0
         )

IF(NOT LIBXMLPP_LIBRARY)
  IF(UNIX)
    IF(RUN_CONFIGURE)
      EXEC_PROGRAM(aclocal ${AUXPACKAGES}/libxml++)
      EXEC_PROGRAM(autoconf ${AUXPACKAGES}/libxml++)
      EXEC_PROGRAM(automake ${AUXPACKAGES}/libxml++ ARGS -a)
      EXEC_PROGRAM(
       "${AUXPACKAGES}/libxml++/configure --prefix=${PROJECT_BINARY_DIR} --with-cxx=$ENV{CXX}"
        ${AUXPACKAGES}/libxml++
      )
      EXEC_PROGRAM(make ${AUXPACKAGES}/libxml++)
      EXEC_PROGRAM("make install" ${AUXPACKAGES}/libxml++)
      FIND_PATH(LIBXMLPP_INCLUDE_DIR libxml++/libxml++.h
                   ${PROJECT_BINARY_DIR}/include/libxml++-1.0)
      FIND_LIBRARY(LIBXMLPP_LIBRARY xml++-0.1
                   ${PROJECT_BINARY_DIR}/lib)
    ENDIF(RUN_CONFIGURE)
  ENDIF(UNIX)
ENDIF(NOT LIBXMLPP_LIBRARY)

INCLUDE_DIRECTORIES(${LIBXMLPP_INCLUDE_DIR})
LINK_LIBRARIES(${LIBXMLPP_LIBRARY})
