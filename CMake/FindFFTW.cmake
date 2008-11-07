#
# this module look for boost (http://www.boost.org) support
# it will define the following values
#
# BOOST_INCLUDE_DIR = where boost/boost.h can be found
#
# May want to define this but seldom required
# BOOST_LIBRARY = where boost library can be found (reserved)
#

IF(FFTW_INCLUDE_DIRS)
  FIND_PATH(FFTW_INCLUDE_DIR fftw3.h  ${FFTW_INCLUDE_DIRS})
  FIND_LIBRARY(FFTW_LIBRARY fftw3 ${FFTW_LIBRARY_DIRS})
ELSE(FFTW_INCLUDE_DIRS)
  SET(TRIAL_PATHS
    $ENV{FFTW_HOME}/include
    /usr/include
    /usr/local/include
    /opt/include
    /usr/apps/include
  )
  FIND_PATH(FFTW_INCLUDE_DIR fftw3.h ${TRIAL_PATHS})

  SET(TRIAL_LIBRARY_PATHS
    $ENV{FFTW_HOME}/lib
    /usr/lib 
    /usr/local/lib
    /opt/lib
    /sw/lib
    )

  FIND_LIBRARY(FFTW_LIBRARY fftw3 ${TRIAL_LIBRARY_PATHS})

ENDIF(FFTW_INCLUDE_DIRS)

IF(FFTW_INCLUDE_DIR)
  SET(FOUND_FFTW 1 CACHE BOOL "Found fftw library")
ELSE(FFTW_INCLUDE_DIR)
  SET(FOUND_FFTW 0 CACHE BOOL "Not Found fftw library")
ENDIF(FFTW_INCLUDE_DIR)

MARK_AS_ADVANCED(
   FFTW_INCLUDE_DIR
   FOUND_FFTW
)
