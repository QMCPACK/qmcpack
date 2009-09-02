# Cmake for cray XT
SET(ENABLE_CATAMOUNT 1)
SET(QMC_BITS 64)
ADD_DEFINITIONS(-DADD_ -DINLINE_ALL=inline -D_CRAYMPI)


IF($ENV{PE_ENV} MATCHES "PGI")
  # pgi environments: no reason to use it !!
  SET(CMAKE_CXX_FLAGS "--restrict -fastsse -Minline=levels:10 --no_exceptions -DBOOST_NO_EXCEPTIONS")

  # link ACML library
  FIND_LIBRARY(BLAS_LIBRARY acml $ENV{ACML_BASE_DIR}/pgi64/lib)
  SET(LAPACK_LIBRARY "")

ELSE($ENV{PE_ENV} MATCHES "PGI")
  IF(CPU_IDENTITY MATCHES "barcelona")
    SET(CMAKE_CXX_FLAGS "-O6 -march=barcelona -msse3 -ftemplate-depth-60 -Drestrict=__restrict__ -fstrict-aliasing -funroll-all-loops   -finline-limit=1000 -ffast-math -Wno-deprecated ")
  ELSE(CPU_IDENTITY MATCHES "barcelona")
    SET(CMAKE_CXX_FLAGS "-O6 -march=opteron -msse3 -ftemplate-depth-60 -Drestrict=__restrict__ -fstrict-aliasing -funroll-all-loops   -finline-limit=1000 -ffast-math -Wno-deprecated ")
  ENDIF(CPU_IDENTITY MATCHES "barcelona")

  SET(HAVE_SSE 1)
  SET(HAVE_SSE2 1)
  SET(HAVE_SSE3 1)
  SET(HAVE_SSSE3 1)
  SET(USE_PREFETCH 1)
  SET(PREFETCH_AHEAD 12)

  FIND_LIBRARY(BLAS_LIBRARY acml $ENV{ACML_BASE_DIR}/gfortran64/lib)
  #SET(LAPACK_LIBRARY "/opt/pgi/6.2.5/linux86-64/6.2/lib/libpgc.a")
  SET(LAPACK_LIBRARY "")
  SET(FORTRAN_LIBS "-lgfortran")

  #openmp is enabled
  SET(CMAKE_TRY_OPENMP_CXX_FLAGS "-fopenmp")
  CHECK_CXX_ACCEPTS_FLAG(${CMAKE_TRY_OPENMP_CXX_FLAGS} GNU_OPENMP_FLAGS)
  IF(GNU_OPENMP_FLAGS)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_TRY_OPENMP_CXX_FLAGS}")
    SET(ENABLE_OPENMP 1)
  ENDIF(GNU_OPENMP_FLAGS)

ENDIF($ENV{PE_ENV} MATCHES "PGI")

IF($ENV{CXX} MATCHES "CC")
  # always link statically
  SET(CMAKE_CXX_LINK_FLAGS " -static")
  SET(HAVE_MPI 1)
  SET(HAVE_OOMPI 1)
  ADD_DEFINITIONS(-DMPICH_SKIP_MPICXX)
ELSE($ENV{CXX} MATCHES "CC")
  SET(ENABLE_CATAMOUNT 0)
ENDIF($ENV{CXX} MATCHES "CC")

## targeting CATAMOUNT
#IF(ENABLE_CATAMOUNT)
#  ADD_DEFINITIONS(-DXT_CATAMOUNT)
#ENDIF(ENABLE_CATAMOUNT)

MARK_AS_ADVANCED(
  LAPACK_LIBRARY 
  BLAS_LIBRARY 
)

##CRAY XT3 CATAMOUNT 
#IF($ENV{CATAMOUNT_DIR} MATCHES "catamount")
#      INCLUDE(${PROJECT_CMAKE}/CrayXTPgi.cmake)
#    ENDIF($ENV{PE_ENV} MATCHES "PGI")
#    IF($ENV{PE_ENV} MATCHES "GNU")
#       INCLUDE(${PROJECT_CMAKE}/GNUCompilers.cmake)
#       IF($ENV{CXX} MATCHES "CC")
#    ENDIF($ENV{PE_ENV} MATCHES "GNU")
#    SET(HAVE_MPI 1)
#    SET(HAVE_OOMPI 1)
#    SET(FOUND_CXXENV 1)
#  ENDIF($ENV{CXX} MATCHES "CC")
#  LINK_LIBRARIES  (/opt/acml/3.6/gnu64/lib/libacml.a;/usr/lib64/libg2c.a)
#ENDIF($ENV{CATAMOUNT_DIR} MATCHES "catamount")

#SET(CMAKE_CXX_FLAGS " --restrict -fast --target=catamount ")
#SET(CMAKE_CXX_FLAGS " --restrict -fast --target=catamount ")
#SET(CMAKE_CXX_FLAGS " --restrict -fastsse -Mipa -Minline=levels:10 --no_exceptions -DBOOST_NO_EXCEPTIONS")
