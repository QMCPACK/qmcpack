
#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER  mpicc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -D_TIMER_ -DUSE_MPI -DMPI_VERSION=3 -D_LINUX_")

set(INTEL_OPTS "-g -malign-double -fomit-frame-pointer -ffast-math -fopenmp -O3 -msse4 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated")
set(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${GNU_OPTS} ${INTEL_OPTS} -std=c++11 ")
set(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${INTEL_OPTS} -std=c99 ")


#--------------------------------------------------------------------------
# below is common for INTEL compilers and MKL library
#--------------------------------------------------------------------------
set(ENABLE_OPENMP 1)
set(HAVE_MPI 1)
set(HAVE_SSE 1)
set(HAVE_SSE2 1)
set(HAVE_SSE3 1)
set(HAVE_SSSE3 1)
set(HAVE_SSE41 1)
set(USE_PREFETCH 1)
set(PREFETCH_AHEAD 10)
set(HAVE_MKL 1)
set(MKL_FOUND 1)
set(HAVE_MKL_VML 1)

set( CMAKE_FIND_ROOT_PATH
 /usr/tce/packages/fftw/fftw-3.3.4-mvapich2-2.2-gcc-4.8-redhat/lib/
 /usr/gapps/qmc/libs/INTEL/boost_1_62_0
 /usr/lib64/ 
 )

# mkl 10.3.x
include_directories(/usr/tce/packages/mkl/mkl-2017.1/include)
set(LAPACK_LIBRARY -L/usr/tce/packages/mkl/mkl-2017.1/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lrt -Wl,-rpath=/usr/tce/packages/mkl/mkl-2017.1/mkl/lib/intel64)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

