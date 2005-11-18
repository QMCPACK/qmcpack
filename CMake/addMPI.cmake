#
# this module look sfor MPI (Message Passing Interface) support
# it will define the following values
#
# MPI_INCLUDE_PATH = where mpi.h can be found
# MPI_LIBRARY      = the library to link against (mpi mpich etc)
#

IF(${CMAKE_SYSTEM_NAME} MATCHES "AIX")

  IF($ENV{MPICXX} MATCHES "mpCC") 
    SET(HAVE_MPI 1)
    SET(HAVE_OOMPI 1)
    MESSAGE(STATUS "Using mpCC to enable MPI")
  ELSE($ENV{MPICXX} MATCHES "mpCC") 
    MESSAGE(FATAL_ERROR "set environment variable CXX to mpCC/mpCC_r for MPI")
  ENDIF($ENV{MPICXX} MATCHES "mpCC") 

ELSE(${CMAKE_SYSTEM_NAME} MATCHES "AIX")

#first, check if vmi is installed: typical linux clusters at NCSA
FIND_PATH(MPI_INCLUDE_PATH mpi.h 
          /usr/local/vmi/mpich/include
)

IF(MPI_INCLUDE_PATH)
  IF(INTEL_COMPILER)
    FIND_LIBRARY(MPI_LIBRARY 
       	         mpich
               	 /usr/local/vmi/mpich/lib/intel 
                 /usr/local/vmi/mpich/lib/icc
                 )
  ELSE(INTEL_COMPILER)
    FIND_LIBRARY(MPI_LIBRARY  mpich
                /usr/local/vmi/mpich/lib/gcc
                )
  ENDIF(INTEL_COMPILER)

  #IF(MPI_LIBRARY)
  #  SET(MPI_LIBRARY ${MPI_LIBRARY} "-lvmi -ldl -lpthread")
  #ENDIF(MPI_LIBRARY)

ENDIF(MPI_INCLUDE_PATH)

MARK_AS_ADVANCED(MPI_INCLUDE_PATH MPI_LIBRARY)

#vmi is not found. Search for the conventional mpi installation
IF(NOT MPI_LIBRARY)
  INCLUDE(${CMAKE_ROOT}/Modules/FindMPI.cmake)
ENDIF(NOT MPI_LIBRARY)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "AIX")

IF(MPI_LIBRARY)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  SET(HAVE_MPI 1)
  SET(HAVE_OOMPI 1)
  SET(FOUND_MPI 1 CACHE BOOL "Found mpi library to link")
ELSE(MPI_LIBRARY)
  IF($ENV{CXX} MATCHES "mpi") 
    SET(HAVE_MPI 1)
    SET(HAVE_OOMPI 1)
    SET(FOUND_MPI 0 CACHE BOOL "No need to link mpi library to link")
  ENDIF($ENV{CXX} MATCHES "mpi") 
ENDIF(MPI_LIBRARY)
