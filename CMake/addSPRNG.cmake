#
# this module look for HDF5 (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# HDF5_INCLUDE_PATH = where hdf5.h can be found
# HDF5_LIBRARY      = the library to link against (hdf5 etc)
#

FIND_LIBRARY(SPRNG_LIBRARY
  sprng
  /usr/local/lib
  /sw/lib
)

FIND_PATH(SPRNG_INCLUDE_PATH
  sprng.h
  /usr/local/include
  /usr/include
  /sw/lib
)

IF(SPRNG_LIBRARY) 
  MESSAGE(STATUS "Found SPRNG library")
  INCLUDE_DIRECTORIES(${SPRNG_INCLUDE_PATH})
ENDIF(SPRNG_LIBRARY) 

LINK_LIBRARIES(${SPRNG_LIBRARY})
ADD_DEFINITIONS(-DUSE_SPRNG)
