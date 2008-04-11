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
#use environment variables: e.g, module at OSC uses MKL
#IF(NOT QMC_BUILD_STATIC)

  IF($ENV{MKL} MATCHES "mkl")
  
    MESSAGE(STATUS "Using intel/mkl library: $ENV{MKL}")
    
    #look for the path where the mkl is installed
    STRING(REGEX MATCHALL "[-][L]([^ ;])+" MKL_PATH_WITH_PREFIX  "$ENV{MKL}")
    STRING(REGEX REPLACE "[-][L]" "" MKL_PATHS ${MKL_PATH_WITH_PREFIX})
  
  ENDIF($ENV{MKL} MATCHES "mkl")
  
  SET(MKL_PATHS ${MKL_PATHS} 
      $ENV{MKL_HOME}/lib/em${QMC_BITS}t
      $ENV{MKL_HOME}/lib/${QMC_BITS}
      $ENV{MKL_HOME}/lib
      /usr/local/intel/mkl60/mkl60/lib/64
      /usr/local/intel/mkl/lib/32 
      /opt/intel/mkl/lib/32
      /opt/intel/mkl/10.0.2.018/lib/em64t
     ) 
  MESSAGE(STATUS "Looking for intel/mkl library in ${MKL_PATHS}")

#ENDIF(NOT QMC_BUILD_STATIC)

IF(NOT LAPACK_LIBRARY_INIT)
  IF(QMC_BITS MATCHES 64)
    FIND_LIBRARY(LAPACK_LIBRARY 
      NAMES mkl_lapack  mkl_ia64
      PATHS ${MKL_PATHS}
      )
  ELSE(QMC_BITS MATCHES 64)
    FIND_LIBRARY(LAPACK_LIBRARY 
      NAMES mkl_lapack mkl_ia32
      PATHS ${MKL_PATHS}
      )
  ENDIF(QMC_BITS MATCHES 64)

  IF(QMC_BUILD_STATIC) 
    FIND_LIBRARY(BLAS_LIBRARY
      NAMES mkl mkl_ipf
      PATHS ${MKL_PATHS}
    )
  ELSE(QMC_BUILD_STATIC) 
    FIND_LIBRARY(BLAS_LIBRARY
      NAMES mkl mkl_itp
      PATHS ${MKL_PATHS}
    )
  ENDIF(QMC_BUILD_STATIC) 

  FIND_LIBRARY(INTEL_GUIDE_LIBRARY
    NAMES guide
    PATHS ${MKL_PATHS}
  )

  IF(LAPACK_LIBRARY MATCHES "mkl")
    MESSAGE(STATUS "Found intel/mkl library")
    SET(LAPACK_LIBRARY_INIT 1 CACHE BOOL "lapack is initialized")
    SET(BLAS_LIBRARY_INIT 1 CACHE BOOL "blas is initialized")
    SET(INTEL_MKL 1 CACHE BOOL "INTEL_MKL is set to 1")
  ENDIF(LAPACK_LIBRARY MATCHES "mkl")

ENDIF(NOT LAPACK_LIBRARY_INIT)

IF(NOT INTEL_MKL)
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    SET(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -framework vecLib")
    SET(LAPACK_LIBRARY_INIT 1 CACHE BOOL "use Mac Framework")
    SET(MAC_VECLIB 1 CACHE BOOL "use Mac Framework")
    SET(LAPACK_LIBRARY "")
    MESSAGE(STATUS "Using Framework on Darwin.")

    ## check goto library: does not work so well
    #FIND_LIBRARY(BLAS_LIBRARY NAMES goto
    #    PATHS 
    #    $ENV{GOTOBLAS_HOME}
    #    /usr/lib
    #    /usr/local/lib
    #    /sw/lib
    #    )
    SET(BLAS_LIBRARY "")
    SET(BLAS_LIBRARY_INIT 1 CACHE BOOL "use Mac Framework")
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
ENDIF(NOT INTEL_MKL)

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

#IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#  SET(LAPACK_LIBRARY $ENV{LAPACK})
#  SET(LAPACK_LIBRARY_INIT 1 CACHE BOOL "lapack is initialized")
#ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

IF(${CMAKE_SYSTEM_NAME} MATCHES "AIX")
  SET(ELIB essl)
  IF(ENABLE_OMP)
    SET(ELIB esslsmp)
  ENDIF(ENABLE_OMP)
   
  IF(NOT LAPACK_LIBRARY_INIT)
    SET(LLIB lapack-SP4_${QMC_BITS} lapack)
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
  FIND_LIBRARY(BLAS_LIBRARY NAMES goto blas blas_gnu
    PATHS 
    $ENV{GOTOBLAS_HOME}
    /usr/apps/math/lapack
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
MARK_AS_ADVANCED(
  LAPACK_LIBRARY 
  BLAS_LIBRARY 
)
