# this module look for lapack/blas and other numerical library support
# it will define the following values
# LAPACK_LIBRARY
# BLAS_LIBRARY
#
# 1) search ENV MKL 
# 2) search MKL in usual paths
# 3) search ENV ATLAS
# 4) search lapack/blas
# 5) give up
# 
SET(MKL_PATHS "")

#use environment variables
IF(INTEL_COMPILER)
  #use environment variables: e.g, module at OSC uses MKL
  IF($ENV{MKL} MATCHES "mkl")

    MESSAGE(STATUS "Using intel/mkl library: $ENV{MKL}")

    #look for the path where the mkl is installed
    STRING(REGEX MATCHALL "[-][L]([^ ;])+" MKL_PATH_WITH_PREFIX  "$ENV{MKL}")
    STRING(REGEX REPLACE "[-][L]" "" MKL_PATHS ${MKL_PATH_WITH_PREFIX})

  ENDIF($ENV{MKL} MATCHES "mkl")

  SET(MKL_PATHS ${MKL_PATHS} 
        $ENV{MKL_HOME}
	/usr/local/intel/mkl60/mkl60/lib/64
	/usr/local/intel/mkl/lib/32 
	/opt/intel/mkl/lib/32
  ) 
  MESSAGE(STATUS "Looking for intel/mkl library in ${MKL_PATHS}")

#  IF($ENV{MKL}) MATCHES "mkl")
#    MESSAGE(STATUS "Using intel/mkl library: $ENV{MKL}")
#    SET(LAPACK_LIBRARY $ENV{MKL})
#    SET(BLAS_LIBRARY " ")
#    SET(LAPACK_LIBRARY_INIT 1)
#  ENDIF($ENV{MKL}) MATCHES "mkl")

  IF(NOT LAPACK_LIBRARY_INIT)
    FIND_LIBRARY(LAPACK_LIBRARY 
      NAMES mkl_lapack 
      PATHS ${MKL_PATHS}
    )
    FIND_LIBRARY(BLAS_LIBRARY
      NAMES mkl mkl_itp
      PATHS ${MKL_PATHS}
    )
    FIND_LIBRARY(BLAS_EXTRA_LIBRARY
      NAMES guide
      PATHS ${MKL_PATHS}
    )
    IF(LAPACK_LIBRARY MATCHES "mkl")
      MESSAGE(STATUS "Found intel/mkl library")
      SET(LAPACK_LIBRARY_INIT 1 CACHE BOOL "lapack is initialized")
      SET(BLAS_LIBRARY_INIT 1 CACHE BOOL "blas is initialized")
    ENDIF(LAPACK_LIBRARY MATCHES "mkl")
  ENDIF(NOT LAPACK_LIBRARY_INIT)
ENDIF(INTEL_COMPILER)

IF($ENV{ATLAS} MATCHES "atlas")
  IF($ENV{ATLAS} MATCHES "lapack") 
    SET(LAPACK_LIBRARY_INIT 1 CACHE BOOL "lapack is initialized with ATLAS env")
  ENDIF($ENV{ATLAS} MATCHES "lapack") 
  SET(BLAS_LIBRARY $ENV{ATLAS})
  SET(BLAS_LIBRARY_INIT 1 CACHE BOOL "blas is initialized with ATLAS env")
ENDIF($ENV{ATLAS} MATCHES "atlas")

IF($ENV{LAPACK} MATCHES "lapack")
  SET(LAPACK_LIBRARY $ENV{LAPACK})
  SET(LAPACK_LIBRARY_INIT 1 CACHE BOOL "lapack is initialized")
ENDIF($ENV{LAPACK} MATCHES "lapack")

IF(${CMAKE_SYSTEM_NAME} MATCHES "AIX")
  SET(ELIB essl)
  IF(ENABLE_OMP)
    SET(ELIB esslsmp)
  ENDIF(ENABLE_OMP)
   
  IF(NOT LAPACK_LIBRARY_INIT)
    SET(LLIB lapack-SP4_32 lapack)
    FIND_LIBRARY(LAPACK_LIBRARY  
                 NAMES ${LLIB}
                 PATHS /usr/apps/math/lapack/LAPACK
                 lib
                )
    FIND_LIBRARY(BLAS_LIBRARY ${ELIB}
                 /usr/lib
                )
  ENDIF(NOT LAPACK_LIBRARY_INIT)

  SET(LAPACK_LIBRARY_INIT 1 CACHE BOOL "lapack is initialized")
  SET(BLAS_LIBRARY_INIT 1 CACHE BOOL "blas with essl is initialized")
  MESSAGE(STATUS "Found lapack/blas on AIX system")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "AIX")

IF(NOT LAPACK_LIBRARY_INIT)
  FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack lapack_gnu
    PATHS /usr/apps/math/lapack
          /usr/lib
          /opt/lib
          /usr/local/lib
          /sw/lib
  )
  IF(LAPACK_LIBRARY)
    MESSAGE(STATUS "Found netlib lapack library")
    SET(LAPACK_LIBRARY_INIT 1 CACHE BOOL "lapack is initialized")
  ENDIF(LAPACK_LIBRARY)
ENDIF(NOT LAPACK_LIBRARY_INIT)

IF(NOT BLAS_LIBRARY_INIT)
  FIND_LIBRARY(BLAS_LIBRARY NAMES blas blas_gnu
    PATHS /usr/apps/math/lapack
          /usr/lib
          /opt/lib
          /usr/local/lib
          /sw/lib
  )
  IF(BLAS_LIBRARY)
    MESSAGE(STATUS "Found netlib blas is found")
    SET(BLAS_LIBRARY_INIT 1 CACHE BOOL "lapack is initialized")
  ENDIF(BLAS_LIBRARY)
ENDIF(NOT BLAS_LIBRARY_INIT)

### BRANDT
### MOVED BECAUSE OF SCOPE PROBLEM
### PROBABLY SHOULD BE FIXED
#SET(BLAS_LIBRARY ${BLAS_LIBRARY} ${BLAS_EXTRA_LIBRARY})
#SET(BLAS_EXTRA_LIBRARY ${BLAS_EXTRA_LIBRARY} ${ATLAS_LIBS})

SET(BLAS_LIBRARY ${BLAS_LIBRARY})
MARK_AS_ADVANCED(
  LAPACK_LIBRARY 
  BLAS_LIBRARY 
)
