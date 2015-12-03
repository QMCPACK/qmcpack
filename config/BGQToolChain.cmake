
set(CMAKE_C_COMPILER mpixlc_r)
set(CMAKE_CXX_COMPILER mpixlcxx_r)

set(AIX_ARCH "qp")
SET(AIX_ARCH_FLAGS "-qarch=${AIX_ARCH} -qsmp=omp -DINLINE_ALL=inline -qthreaded -qstrict -qhot=level=1 -qtune=qp -qsimd=auto  -DHAVE_MASS -DHAVE_MASSV -DBGQPX -DSPLINEFLOAT -DUSE_REAL_STRUCT_FACTOR")

SET(AIX_CXX_COMMON_FLAGS "-g -O3 -qinline=auto:level=10")
SET(AIX_OPT_FLAGS "-qmaxmem=-1")

SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS}")


SET(BGQ 1)
SET(QMC_CUDA 0)
SET(QMC_COMPLEX 0)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
set(HAVE_CUDA 0)
SET(QMC_BUILD_STATIC 1)
SET(HAVE_LIBESSL 1)
SET(HAVE_EINSPLINE 1)
SET(HAVE_EINSPLINE_EXT 0)
SET(HAVE_ADIOS 0)
SET(BUILD_QMCTOOLS 1)
SET(BUILD_SANDBOX 1)




SET(Boost_INCLUDE_DIR /home/projects/qmcpack/boost_1_45_0)
SET(CMAKE_FIND_ROOT_PATH
     /home/projects/qmcpack/boost_1_45_0
     /home/projects/qmcpack/LIBXML2-2.9
     /soft/libraries/hdf5/current/cnk-xl/current
     /soft/libraries/alcf/current/xl/FFTW3
     /soft/libraries/alcf/current/xl/ZLIB
)


SET(LAPACK_LIBRARY /soft/libraries/alcf/current/xl/LAPACK/lib/liblapack.a)
SET(BLAS_LIBRARY /soft/libraries/essl/current/essl/5.1/lib64/libesslsmpbg.a)
SET(FORTRAN_LIBRARIES
$ENV{IBM_FCMP_DIR}/bglib64/libxlf90_r.a
$ENV{IBM_FCMP_DIR}/bglib64/libxlopt.a
)
LINK_LIBRARIES(
$ENV{IBM_MAIN_DIR}/xlmass/bg/7.3/bglib64/libmass.a 
$ENV{IBM_MAIN_DIR}/xlmass/bg/7.3/bglib64/libmassv.a 
#/soft/perftools/hpctw/libmpihpm_smp.a
#/bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a
#/bgsys/drivers/ppcfloor/spi/lib/libSPI_upci_cnk.a
#-pg
)

