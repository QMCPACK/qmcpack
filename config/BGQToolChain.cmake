# the name of the target operating system
#SET(CMAKE_SYSTEM_NAME BlueGeneP)
SET(BGP 1 CACHE BOOL "On BlueGeneP")
SET(Linux 0)

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_LIBESSL 1)

# set the compiler
set(CMAKE_C_COMPILER /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlc_r)
set(CMAKE_CXX_COMPILER /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx_r)

set(AIX_ARCH "qp")
SET(AIX_ARCH_FLAGS "-qarch=${AIX_ARCH} -qsmp=omp -qthreaded -qstrict")
SET(AIX_CXX_COMMON_FLAGS "-g -O3 -qinline")
SET(AIX_OPT_FLAGS "-qmaxmem=-1")

SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS}")

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(Boost_INCLUDE_DIR /home/xanaromero/boost_1_45_0)
set(CMAKE_FIND_ROOT_PATH
    /home/xanaromero/boost_1_45_0
    /home/xanaromero/einspline-0.9.2
    /home/xanaromero/libxml2-2.7.2
    /soft/apps/hdf5-1.8.7
    /soft/apps/fftw3/fftw-3.3-xl
)

SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)

set(LAPACK_LIBRARY /soft/apps/netlib/liblapack_BGQ_XL.a)
set(BLAS_LIBRARY /soft/apps/essl/prerelease/libesslbg.a)
set(FORTRAN_LIBRARIES 
/opt/ibmcmp/xlf/bg/14.1/bglib64/libxlf90_r.a
/opt/ibmcmp/xlf/bg/14.1/bglib64/libxlfmath.a
/opt/ibmcmp/xlf/bg/14.1/bglib64/libxlopt.a
)

SET(BUILD_SANDBOX 1)
