# common variables for power5+
# disable some inline keywords
# attach the flags which are used by external libraries 
ADD_DEFINITIONS(-DINLINE_ALL= -DIBM)

SET(HAVE_ESSL 1)
SET(BLAS_LIBRARY -lessl -lmassv -lmass)
SET(LAPACK_LIBRARY -llapack)

IF(${CMAKE_CXX_COMPILER} MATCHES "mp")
  SET(HAVE_MPI 1)
  SET(HAVE_OOMPI 1)
ENDIF(${CMAKE_CXX_COMPILER} MATCHES "mp")

IF(${CMAKE_C_COMPILER} MATCHES "_r")
  SET(ENABLE_OPENMP 1)
ENDIF(${CMAKE_C_COMPILER} MATCHES "_r")
SET(FORTRAN_LIBS  -L/usr/lpp/xlf/lib -lxlf90 -lxlf -lxlopt)

MARK_AS_ADVANCED(
  LAPACK_LIBRARY 
  BLAS_LIBRARY 
)
