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

SET(MKL_PATHS "/usr/local/lib /usr/lib")

#IF(NOT QMC_BUILD_STATIC)
IF($ENV{MKL} MATCHES "mkl")
  MESSAGE(STATUS "Using intel/mkl library: $ENV{MKL}")
  #look for the path where the mkl is installed
  STRING(REGEX MATCHALL "[-][L]([^ ;])+" MKL_PATH_WITH_PREFIX  "$ENV{MKL}")
  STRING(REGEX REPLACE "[-][L]" "" MKL_PATHS ${MKL_PATH_WITH_PREFIX})
ENDIF($ENV{MKL} MATCHES "mkl")

SET(MKL_PATHS 
  $ENV{MKL_HOME}/lib/em${QMC_BITS}t
  $ENV{MKL_HOME}/lib/${QMC_BITS}
  $ENV{MKL_HOME}/lib
  ${MKL_PATHS} 
  ) 
MESSAGE(STATUS "Looking for intel/mkl library in ${MKL_PATHS}")

STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" MKL_VERSION "$ENV{MKL_HOME}")
SET(LINK_LAPACK_DEFAULT 1)
IF(${MKL_VERSION} MATCHES "10\\.0\\.[0-2]")
  SET(LINK_LAPACK_DEFAULT 0)
ENDIF(${MKL_VERSION} MATCHES "10\\.0\\.[0-2]")

IF(NOT LAPACK_LIBRARY_INIT)
  IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    IF(LINK_LAPACK_DEFAULT)
      FIND_LIBRARY(LAPACK_LIBRARY NAMES mkl_lapack PATHS ${MKL_PATHS})
      FIND_LIBRARY(BLAS_LIBRARY NAMES mkl PATHS ${MKL_PATHS})
    ELSE(LINK_LAPACK_DEFAULT)
      FIND_LIBRARY(LAPACK_LIBRARY STATIC NAMES mkl_lapack PATHS ${MKL_PATHS})
      FIND_LIBRARY(BLAS_LIBRARY  NAMES mkl_em64t PATHS ${MKL_PATHS})
      MESSAGE("-- mkl 10.0.[0-2] warning for EM64T")
      MESSAGE("-- Replace libmkl_lapack.so in CMakeCache.txt by libmkl_lapack.a")
    ENDIF(LINK_LAPACK_DEFAULT)
  ELSE(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    FIND_LIBRARY(LAPACK_LIBRARY 
      NAMES mkl_lapack 
      PATHS ${MKL_PATHS}
      )
    FIND_LIBRARY(BLAS_LIBRARY
      NAMES mkl
      PATHS ${MKL_PATHS}
    )
  ENDIF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")

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

IF(USE_SCALAPACK)
  SET(PNPATHS 
    ${MKL_PATHS}
    ${BLACS_HOME}/lib
    ${SCALAPACK_HOME}/lib
    /usr/lib
    /opt/lib
    /usr/local/lib
    /sw/lib
    )

  IF(INTEL_MKL)
    FIND_LIBRARY(BLACSLIB mkl_blacs_${PLAT}_lp${QMC_BITS} PATHS  ${PNPATHS})
    FIND_LIBRARY(SCALAPACKLIB mkl_scalapack PATHS  ${PNPATHS})
  ENDIF(INTEL_MKL)

  IF(NOT SCALAPACKLIB)
    FIND_LIBRARY(BLACSLIB blacs_MPI-${PLAT}-{BLACSDBGLVL} PATHS  ${PNPATHS})
    FIND_LIBRARY(BLACSCINIT blacsCinit_MPI-${PLAT}-{BLACSDBGLVL} PATHS  ${PNPATHS})
    FIND_LIBRARY(SCALAPACKLIB scalapack PATHS  ${PNPATHS})
  ENDIF(NOT SCALAPACKLIB)

  IF(BLACSLIB AND SCALAPACKLIB)
    SET(FOUND_SCALAPACK 1 CACHE BOOL "Found scalapack library")
  ELSE(BLACSLIB AND SCALAPACKLIB)
    SET(FOUND_SCALAPACK 0 CACHE BOOL "Mising scalapack library")
  ENDIF(BLACSLIB AND SCALAPACKLIB)

  MARK_AS_ADVANCED(
    BLACSCINIT
    BLACSLIB
    SCALAPACKLIB
    FOUND_SCALAPACK
    )
ENDIF(USE_SCALAPACK)

