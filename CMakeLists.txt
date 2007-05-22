PROJECT(qmcpack)

######################################################################
# Version 
######################################################################
SET(QMCPLUSPLUS_VERSION_MAJOR 0)
SET(QMCPLUSPLUS_VERSION_MINOR 2)
SET(QMCPLUSPLUS_VERSION_PATCH 0)
SET(QMCPLUSPLUS_RELEASE 0)
SET(QMCPLUSPLUS_VERSION "${QMCPLUSPLUS_VERSION_MAJOR}.${QMCPLUSPLUS_VERSION_MINOR}.$(QMCPLUSPLUS_VERSION_PATH}")

######################################################################
#define the paths for library and executable for out-of-soruce compilation
######################################################################
SET (LIBRARY_OUTPUT_PATH ${qmcpack_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${qmcpack_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

######################################################################
# Set the compiler-time parameters
# OHMMS_DIM =  dimension of the problem
# OHMMS_INDEXTYPE = type of index
# OHMMS_PRECISION  = base precision, float, double etc
# QMC_COMPLEX = true if using complex wavefunctions
# QMC_MPI =  enable MPI 
# QMC_OMP = enable OMP 
# QMC_BITS = 32/64 bit
######################################################################
SET(OHMMS_DIM 3)
SET(OHMMS_INDEXTYPE int)
SET(OHMMS_PRECISION double)
SET(QMC_COMPLEX 0)
SET(QMC_MPI 1)
SET(QMC_OMP 0)
SET(QMC_BITS 32)

######################################################################
# Standard test
######################################################################
INCLUDE(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
CHECK_FUNCTION_EXISTS(sincos HAVE_SINCOS)

######################################################################
# Build options
# QMC_BUILD_STATIC build static/dynamic  libraries
# BUILD_QMCTOOLS   build utility programs
# BUILD_SANDBOX    build test programs
# MPIP_PROFILE     profile mpi performance
######################################################################
SET(QMC_BUILD_STATIC 0)
SET(BUILD_QMCTOOLS 0)
SET(BUILD_SANDBOX 0)
SET(MPIP_PROFILE 0)

######################################################################
# Performance-related macros
# QMC_SK_USE_RECURSIVE enable/disable recursive evalaution of SK
######################################################################
SET(QMC_SK_USE_RECURSIVE 0)

######################################################################
# Overwrite with environments if they differ from the set values
######################################################################
IF(NOT $ENV{QMC_MPI} MATCHES QMC_MPI)
  SET(QMC_MPI $ENV{QMC_MPI}) 
ENDIF(NOT $ENV{QMC_MPI} MATCHES QMC_MPI)
IF(NOT $ENV{QMC_OMP} MATCHES QMC_OMP)
  SET(QMC_OMP $ENV{QMC_OMP}) 
ENDIF(NOT $ENV{QMC_OMP} MATCHES QMC_OMP)
IF(NOT $ENV{QMC_BITS} MATCHES QMC_BITS)
  SET(QMC_BITS $ENV{QMC_BITS})
ENDIF(NOT $ENV{QMC_BITS} MATCHES QMC_BITS)
IF(NOT $ENV{QMC_COMPLEX} MATCHES QMC_COMPLEX)
  SET(QMC_COMPLEX $ENV{QMC_COMPLEX})
ENDIF(NOT $ENV{QMC_COMPLEX} MATCHES QMC_COMPLEX)

######################################################################
# FIXED PARAMETERS for test and legacy reasons
# DO NOT TOUCH THESE
# QMC_ORTHO  = 0 for general unit cell, 1 for orthogonal cells
# QMC_FASTWALKER  = true, when operation is faster for each walker
######################################################################
SET(QMC_FASTWALKER 1)
SET(OHMMS_ORTHO 0)

######################################################################
#list of options to determine libraries. First, try to use libraries available
######################################################################
SET(INTEL_COMPILER 0)
SET(ENABLE_OPENMP 0)
SET(HAVE_MPI 0)
SET(HAVE_OOMPI 0)
SET(HAVE_LIBXML2 1)
SET(HAVE_LIBFFTW 0)
SET(HAVE_LIBXMLPP 0)
SET(HAVE_LIBHDF5 0)
SET(HAVE_LIBSPRNG 0)
SET(HAVE_LIBBLITZ 0)
SET(HAVE_LIBGSL 0)
SET(INTEL_COMPILER 0)
SET(MAC_VECLIB 0)

######################################################
# Directory where customize cmake files reside
######################################################
SET (PROJECT_CMAKE ${qmcpack_SOURCE_DIR}/CMake)

######################################################
#COMPILER choose one of the cmake files to customize the compiler options
#If nothing is chosen, default settings by cmake will  be used.
#If the automatic detection does not work, comment out everything 
#upto COMPILER MANUAL SELECTION and use one of the customized cmake file.
######################################################
SET(FOUND_CXXENV 0)
IF($ENV{CXX} MATCHES "icpc")
  INCLUDE(${PROJECT_CMAKE}/Intel8.cmake)
  SET(FOUND_CXXENV 1)
ENDIF($ENV{CXX} MATCHES "icpc")
IF($ENV{CXX} MATCHES "cmpic")
  INCLUDE(${PROJECT_CMAKE}/Intel8.cmake)
  SET(FOUND_CXXENV 1)
ENDIF($ENV{CXX} MATCHES "cmpic")
IF($ENV{CXX} MATCHES "xlC")
  INCLUDE(${PROJECT_CMAKE}/IBMCompilers.cmake)
  SET(FOUND_CXXENV 1)
ENDIF($ENV{CXX} MATCHES "xlC")

IF(NOT FOUND_CXXENV)
  INCLUDE(${PROJECT_CMAKE}/GNUCompilers.cmake)
ENDIF(NOT FOUND_CXXENV)

#COMPIELR MANUAL SELECTION
#INCLUDE(${PROJECT_CMAKE}/GNUCompilers.cmake)
#INCLUDE(${PROJECT_CMAKE}/Intel8.cmake)
#INCLUDE(${PROJECT_CMAKE}/IBMCompilers.cmake)
#INCLUDE(${PROJECT_CMAKE}/KCCCompilers.cmake)
#INCLUDE(${PROJECT_CMAKE}/Intel7.cmake)

######################################################
# DO NOT MODIFY BELOW
######################################################
INCLUDE(${CMAKE_ROOT}/Modules/FindThreads.cmake)
INCLUDE(${PROJECT_CMAKE}/FindLibxml2.cmake)
INCLUDE(${PROJECT_CMAKE}/FindLapack.cmake)
INCLUDE(${PROJECT_CMAKE}/FindHDF5.cmake)
INCLUDE(${PROJECT_CMAKE}/FindBoost.cmake)
#INCLUDE(${PROJECT_CMAKE}/FindFFTW.cmake)

IF(QMC_MPI)
  IF(MPIP_PROFILE)
    LINK_LIBRARIES("/usr/apps/tools/mpip/lib/cmpi/libmpiP.a")
    LINK_LIBRARIES("/usr/lib/libbfd.a")
    LINK_LIBRARIES("/usr/lib/libiberty.a")
  ENDIF(MPIP_PROFILE)

  INCLUDE(${PROJECT_SOURCE_DIR}/CMake/addMPI.cmake)
  IF(FOUND_MPI)
    LINK_LIBRARIES(${MPI_LIBRARY})
  ENDIF(FOUND_MPI)
ENDIF(QMC_MPI)

IF(FOUND_LIBXML2)
  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
  LINK_LIBRARIES(${LIBXML2_LIBRARY})
  SET(HAVE_LIBXML2 1)
ENDIF(FOUND_LIBXML2)

IF(FOUND_BOOST)
  SET(HAVE_LIBBOOST 1)
  INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})
ENDIF(FOUND_BOOST)

#IF(FOUND_FFTW)
#  SET(HAVE_LIBFFTW 1)
#  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
#  LINK_LIBRARIES(${FFTW_LIBRARY})
#ENDIF(FOUND_FFTW)

IF(FOUND_HDF5)
  SET(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  LINK_LIBRARIES(${HDF5_LIBRARY})
  #zlib
  NCLUDE(${PROJECT_CMAKE}/FindZLIB.cmake)
  IF(ZLIB_FOUND)
    LINK_LIBRARIES(${ZLIB_LIBRARY})
  ENDIF(ZLIB_FOUND)
  #szlib
#  INCLUDE(${PROJECT_CMAKE}/FindSZLIB.cmake)
#  IF(SZLIB_FOUND)
#    LINK_LIBRARIES(${SZLIB_LIBRARY})
#  ENDIF(SZLIB_FOUND)
ENDIF(FOUND_HDF5)

INCLUDE_DIRECTORIES(
  ${qmcpack_SOURCE_DIR}/src
  ${qmcpack_BINARY_DIR}/src
)
ADD_DEFINITIONS(-DHAVE_CONFIG_H -DAPPNAMESPACE=qmcplusplus)

SUBDIRS(src)
