#
# this module look for GSL (http://sources.redhat.com/gsl/) support
# it will define the following values
#
# GSL_INCLUDE_PATH = where gsl/gsl_multimin.h can be found
# GSL_LIBRARY      = the library to link against (gsl etc)
#

FIND_LIBRARY(GSL_LIBRARY
  gsl  
  /u/ac/esler/lib/gsl/lib
  /usr/lib
  /opt/lib
  /usr/local/lib
  /sw/lib
)
FIND_LIBRARY(GSLCBLAS_LIBRARY
  gslcblas
  /u/ac/esler/lib/gsl/lib
  /usr/lib
  /opt/lib
  /usr/local/lib
  /sw/lib
)
FIND_PATH(GSL_INCLUDE_PATH
  gsl/gsl_vector.h
  /u/ac/esler/lib/gsl/include
  /usr/include
  /opt/include
  /usr/local/include
  /sw/include
)

IF(GSL_LIBRARY) 
  MESSAGE(STATUS "Found GSL")
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_PATH})
  SET(GSL_LIBRARY ${GSL_LIBRARY})
ELSE(GSL_LIBRARY) 
  INCLUDE_DIRECTORIES(${AUXPACKAGES}/Utilities/gsl ${AUXPACKAGES}/Utilities/gsl)
  SUBDIRS(Utilities)
ENDIF(GSL_LIBRARY) 

IF(GSLCBLAS_LIBRARY) 
  MESSAGE(STATUS "Found GSLCBLAS")
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_PATH})
  SET(GSLCBLAS_LIBRARY ${GSLCBLAS_LIBRARY})
ELSE(GSLCBLAS_LIBRARY) 
  INCLUDE_DIRECTORIES(${AUXPACKAGES}/Utilities/gslcblas ${AUXPACKAGES}/Utilities/gslcblas)
  SUBDIRS(Utilities)
ENDIF(GSLCBLAS_LIBRARY) 

LINK_LIBRARIES(${GSL_LIBRARY})
LINK_LIBRARIES(${GSLCBLAS_LIBRARY})
