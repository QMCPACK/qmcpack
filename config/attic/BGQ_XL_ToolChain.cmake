
set(CMAKE_C_COMPILER mpixlc_r)
set(CMAKE_CXX_COMPILER mpixlcxx_r)

# set the linker if extra wrapper is needed for profiling tools like hpctoolkit.
#set(CMAKE_CXX_LINKER "hpclink mpixlcxx_r")

set(AIX_ARCH "qp")
SET(AIX_ARCH_FLAGS "-qarch=${AIX_ARCH} -qsmp=omp -qthreaded -qstrict -qhot=level=1 -qtune=qp -qsimd=auto -DHAVE_MASS -DHAVE_MASSV -DBGQPX -DSPLINEFLOAT -D__forceinline=inline")

SET(AIX_SUPPRESSED_WARN "-qsuppress=1540-0724:1540-0198:1540-1401:1540-0822:1540-0095:1540-1101:1500-010:1500-029")
SET(AIX_CXX_COMMON_FLAGS "-g -O3 -qinline=auto:level=10 ${AIX_SUPPRESSED_WARN}")
SET(AIX_OPT_FLAGS "-qmaxmem=-1")

SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS} ${AIX_CXX_COMMON_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS}")

SET(QMC_CUDA 0)
#SET(QMC_COMPLEX 0)
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

SET(MPIEXEC "sh")
SET(MPIEXEC_NUMPROC_FLAG "${qmcpack_SOURCE_DIR}/utils/bgrunjobhelper.sh")
SET(QE_BIN /soft/applications/quantum_espresso/5.3.0-bgq-omp/bin)

SET(BOOST_ROOT /home/projects/qmcpack/boost_1_45_0)

SET(CMAKE_FIND_ROOT_PATH
     /home/projects/qmcpack/libXML2-2.9.1
     /soft/libraries/hdf5/1.8.14/cnk-xl/current
     /soft/libraries/alcf/current/xl/FFTW3
     /soft/libraries/alcf/current/xl/ZLIB
)

SET(LAPACK_LIBRARY /soft/libraries/alcf/current/xl/LAPACK/lib/liblapack.a)
SET(BLAS_LIBRARY /soft/libraries/essl/current/essl/5.1/lib64/libesslsmpbg.a)
SET(FORTRAN_LIBRARIES
$ENV{IBM_MAIN_DIR}/xlmass/bg/7.3/bglib64/libmass.a 
$ENV{IBM_MAIN_DIR}/xlmass/bg/7.3/bglib64/libmassv.a 
$ENV{IBM_FCMP_DIR}/bglib64/libxlf90_r.a
$ENV{IBM_FCMP_DIR}/bglib64/libxlopt.a
)

#LINK_LIBRARIES(
#/soft/perftools/hpctw/libmpihpm_smp.a
#/bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a
#/bgsys/drivers/ppcfloor/spi/lib/libSPI_upci_cnk.a
#-pg
#)

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)

