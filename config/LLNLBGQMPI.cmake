# the name of the target operating system
#SET(CMAKE_SYSTEM_NAME BlueGeneP)
SET(BGP 1 CACHE BOOL "On BlueGeneQ")
SET(Linux 0)

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_LIBESSL 1)
SET(HAVE_EINSPLINE 1)
SET(HAVE_EINSPLINE_EXT 0)

# set the compiler
set(CMAKE_C_COMPILER mpixlc_r)
set(CMAKE_CXX_COMPILER mpixlcxx_r)

set(AIX_ARCH "qp")
SET(AIX_ARCH_FLAGS "-qarch=${AIX_ARCH} -qsmp=omp -DINLINE_ALL=inline -qthreaded -qstrict -qhot=level=1 -qtune=qp -qsimd=auto  -DHAVE_MASS -DHAVE_MASSV -DBGQPX -DSPLINEFLOAT -DUSE_REAL_STRUCT_FACTOR")

SET(AIX_CXX_COMMON_FLAGS "-g -O3 -qinline=auto:level=10")
SET(AIX_OPT_FLAGS "-qmaxmem=-1")

SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS}")

set(Boost_INCLUDE_DIR /usr/global/tools/adept/boost-1.46.1/bgqos_0/include)
set(Boost_LIBRARY_DIRS /usr/global/tools/adept/boost-1.46.1/bgqos_0/lib/boost-1.46.1)
set(CMAKE_FIND_ROOT_PATH
    /usr/local/tools/hdf5/hdf5-1.8.5/serial/
    /usr/local/tools/fftw-3.3.1/
    /usr/gapps/qmc/libs/BGQ/libxml2-2.7.4/
)


SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)

set(LAPACK_LIBRARY /usr/local/tools/lapack/lib/liblapack.a) 
set(BLAS_LIBRARY /usr/local/tools/essl/5.1/lib/libesslsmpbg.a)
set(FORTRAN_LIBRARIES 
/opt/ibmcmp/xlf/bg/14.1/bglib64/libxlf90_r.a
#/opt/ibmcmp/xlf/bg/14.1/bglib64/libxlfmath.a
/opt/ibmcmp/xlf/bg/14.1/bglib64/libxlopt.a
)

link_libraries(
/opt/ibmcmp/xlmass/bg/7.3/bglib64/libmass.a 
/opt/ibmcmp/xlmass/bg/7.3/bglib64/libmassv.a 
#/soft/perftools/hpctw/libmpihpm_smp.a
#/bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a
#/bgsys/drivers/ppcfloor/spi/lib/libSPI_upci_cnk.a
#-pg
)
SET(BUILD_QMCTOOLS 1)
SET(BUILD_SANDBOX 1)
