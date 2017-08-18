set(CMAKE_C_COMPILER mpiclang) 
set(CMAKE_CXX_COMPILER mpiclang++11)

set(GNU_OPTS "-O3 -g -ffast-math -fopenmp -fstrict-aliasing -Wno-deprecated -Wno-unused-value -Wno-type-safety -Wno-undefined-var-template")
set(GNU_FLAGS "-Drestrict=__restrict__ -DADD_ -DHAVE_MASS -DHAVE_MASSV -DSPLINEFLOAT -DBGQPX -D__forceinline=inline")
set(CMAKE_CXX_FLAGS "${GNU_FLAGS} ${GNU_OPTS} -ftemplate-depth-60")
set(CMAKE_C_FLAGS "${GNU_FLAGS} ${GNU_OPTS} -std=c99" )
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")

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
#SET(BUILD_SANDBOX 1)

SET(MPIEXEC "sh")
SET(MPIEXEC_NUMPROC_FLAG "${qmcpack_SOURCE_DIR}/utils/bgrunjobhelper.sh")
SET(QE_BIN /soft/applications/quantum_espresso/5.3.0-bgq-omp/bin)

SET(BOOST_ROOT /soft/libraries/boost/1.62.0/cnk-bgclang++11/current/)

SET(CMAKE_FIND_ROOT_PATH
     /home/projects/qmcpack/libXML2-2.9.1
     /soft/libraries/hdf5/current/cnk-gcc/current
     /soft/libraries/alcf/current/gcc/FFTW3
     /soft/libraries/alcf/current/gcc/ZLIB
)

include_directories($ENV{IBM_MAIN_DIR}/xlmass/bg/7.3/include)

SET(LAPACK_LIBRARY /soft/libraries/alcf/current/xl/LAPACK/lib/liblapack.a)
SET(BLAS_LIBRARY /soft/libraries/essl/current/essl/5.1/lib64/libesslbg.a)
SET(FORTRAN_LIBRARIES
$ENV{IBM_MAIN_DIR}/xlmass/bg/7.3/bglib64/libmass.a 
$ENV{IBM_MAIN_DIR}/xlmass/bg/7.3/bglib64/libmassv.a 
$ENV{IBM_FCMP_DIR}/bglib64/libxlf90_r.a
$ENV{IBM_FCMP_DIR}/bglib64/libxlopt.a
$ENV{IBM_FCMP_DIR}/bglib64/libxl.a
/soft/compilers/bgclang/xlsmp-nonconflicting/ibmcmp-feb2015/libxlsmp.a
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

