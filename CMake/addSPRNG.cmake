#
# this module look for HDF5 (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# HDF5_INCLUDE_PATH = where hdf5.h can be found
# HDF5_LIBRARY      = the library to link against (hdf5 etc)
#

SET(TRIAL_LIBRARY_PATHS
  /usr/local/lib
  /sw/lib
)
SET(TRIAL_INCLUDE_PATHS
  /usr/local/include
  /sw/include
)

IF($ENV{SPRNG_HOME} MATCHES "sprng")
  SET(TRIAL_LIBRARY_PATHS
     $ENV{SPRNG_HOME}/lib ${TRIAL_LIBRARY_PATHS} )
  SET(TRIAL_INCLUDE_PATHS
     $ENV{SPRNG_HOME}/include ${TRIAL_INCLUDE_PATHS} )
ENDIF($ENV{SPRNG_HOME} MATCHES "sprng")

FIND_LIBRARY(SPRNG_LIBRARY
  sprng
  ${TRIAL_LIBRARY_PATHS}
)

FIND_PATH(SPRNG_INCLUDE_PATH
  sprng.h
  ${TRIAL_INCLUDE_PATHS} 
)

IF(SPRNG_LIBRARY) 
  MESSAGE(STATUS "Found SPRNG library")
  INCLUDE_DIRECTORIES(${SPRNG_INCLUDE_PATH})
ENDIF(SPRNG_LIBRARY) 

LINK_LIBRARIES(${SPRNG_LIBRARY})
ADD_DEFINITIONS(-DUSE_SPRNG)
