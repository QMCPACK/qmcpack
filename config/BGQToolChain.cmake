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

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(Boost_INCLUDE_DIR /home/projects/qmcpack/boost_1_45_0)
SET(CMAKE_FIND_ROOT_PATH
     /home/projects/qmcpack/boost_1_45_0
     /home/projects/qmcpack/LIBXML2-2.9
     /soft/libraries/unsupported/hdf5-1.8.8
     /soft/libraries/alcf/current/xl/FFTW3
)

SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)

set(LAPACK_LIBRARY /soft/libraries/alcf/current/xl/LAPACK/lib/liblapack.a)
set(BLAS_LIBRARY /soft/libraries/essl/5.1.1-0/essl/5.1/lib64/libesslsmpbg.a)
SET(FORTRAN_LIBRARIES
/soft/compilers/ibmcmp-feb2013/xlf/bg/14.1/bglib64/libxlf90_r.a
/soft/compilers/ibmcmp-feb2013/xlf/bg/14.1/bglib64/libxlopt.a
)
link_libraries(
/soft/compilers/ibmcmp-feb2013/xlmass/bg/7.3/bglib64/libmass.a 
/soft/compilers/ibmcmp-feb2013/xlmass/bg/7.3/bglib64/libmassv.a 
#/soft/perftools/hpctw/libmpihpm_smp.a
#/bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a
#/bgsys/drivers/ppcfloor/spi/lib/libSPI_upci_cnk.a
#-pg
)
SET(BUILD_QMCTOOLS 1)
SET(BUILD_SANDBOX 1)
