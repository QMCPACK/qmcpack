# the name of the target operating system
#SET(CMAKE_SYSTEM_NAME BlueGeneP)
SET(BGP 1 CACHE BOOL "On BlueGeneP")
SET(Linux 0)

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_LIBESSL 1)

# set the compiler
#set(CMAKE_C_COMPILER  /opt/ibmcmp/vacpp/bg/9.0/bin/bgxlc_r)
#set(CMAKE_CXX_COMPILER  /opt/ibmcmp/vacpp/bg/9.0/bin/bgxlC_r)
set(CMAKE_C_COMPILER  /soft/apps/current/gcc-4.3.2/comm/default/bin/mpicc)             
set(CMAKE_CXX_COMPILER  /soft/apps/current/gcc-4.3.2/comm/default/bin/mpicxx )           

SET(CMAKE_CXX_FLAGS "-g -O2 -Drestrict=__restrict__  -Wno-deprecated  -fopenmp")
SET(CMAKE_C_FLAGS "-O2 -g -Drestrict=__restrict__  -std=gnu99 -fomit-frame-pointer -fopenmp ")

## use this to search only these directores not standard directories.
SET(CMAKE_FIND_ROOT_PATH  
/home/projects/qmcpack/libxml_gcc
/home/projects/qmcpack/boost-1.45_0
/home/projects/qmcpack/einspline_gcc
/soft/apps/hdf5-1.8.0
/soft/apps/fftw-3.1.2-double
)
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)

set(LAPACK_LIBRARY /home/projects/qmcpack/liblapack_gcc.a)
set(BLAS_LIBRARY /soft/apps/ESSL-4.4.1-1/lib/libesslbg.a)
SET(FORTRAN_LIBRARIES 
/soft/apps/ibmcmp-aug2011/xlf/bg/11.1/bglib/libxlf90_r.a
/soft/apps/ibmcmp-aug2011/xlf/bg/11.1/bglib/libxlfmath.a
/soft/apps/ibmcmp-aug2011/xlf/bg/11.1/bglib/libxlopt.a
/soft/apps/ibmcmp-aug2011/xlf/bg/11.1/bglib/libxl.a
)

