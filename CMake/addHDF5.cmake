#
# this module look for HDF5 (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# HDF5_INCLUDE_PATH = where hdf5.h can be found
# HDF5_LIBRARY      = the library to link against (hdf5 etc)
#

FIND_LIBRARY(HDF5_LIBRARY
  hdf5
  /home/common/lib/hdf5/lib
  /u/ac/esler/lib/hdf5/lib
  /usr/lib
  /opt/lib
  /usr/local/lib
  /usr/apps/lib
  /sw/lib
)
FIND_PATH(HDF5_INCLUDE_PATH
  H5Ipublic.h
  /home/common/lib/hdf5/include
  /u/ac/esler/lib/hdf5/include
  /usr/lib
  /opt/lib
  /usr/local/include
  /usr/apps/include
  /sw/lib
)

IF(HDF5_LIBRARY) 
  MESSAGE(STATUS "Found HDF5")
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_PATH})
  SET(HDF5_LIBRARY ${HDF5_LIBRARY})
ELSE(HDF5_LIBRARY) 
  INCLUDE_DIRECTORIES(${AUXPACKAGES}/Utilities/hdf5 ${AUXPACKAGES}/Utilities/hdf5)
  SUBDIRS(Utilities)
ENDIF(HDF5_LIBRARY) 

LINK_LIBRARIES(${HDF5_LIBRARY})
