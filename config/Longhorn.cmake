#--------------------------------------------------------------------------
# tool chain for abe Wed Aug 12 2009
#--------------------------------------------------------------------------
SET(CMAKE_SYSTEM_PROCESSOR "ES")

#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER /opt/apps/intel11_1/mvapich2/1.4/bin/mpicxx)
set(CMAKE_C_COMPILER  /opt/apps/intel/11.1/bin/intel64/icc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(INTEL_OPTS "-g  -restrict -unroll  -O3 -ip -xSSE4.2 -openmp -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS} -std=c99")
set(MKL_HOME "/usr/local/intel/mkl/10.1.2.024" CACHE STRING "MKL HOME")

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(EINSPLINE_HOME "/u/ac/jnkim/svnwork/einspline" CACHE STRING "Einspline source directory")
set(CMAKE_FIND_ROOT_PATH
    /home/00504/tg457645/build/boost_1_42_0
    /opt/apps/intel11_1/mvapich2_1_4/hdf5/1.8.3
    /home/00504/tg457645/
    )

#--------------------------------------------------------------------------
# below is common for INTEL compilers and MKL library
#--------------------------------------------------------------------------
set(ENABLE_OPENMP 1)
set(HAVE_MPI 1)
set(HAVE_SSE 1)
set(HAVE_SSE2 1)
set(HAVE_SSE3 1)
set(HAVE_SSSE3 1)
set(USE_PREFETCH 1)
set(PREFETCH_AHEAD 10)
set(HAVE_MKL 0)
set(HAVE_MKL_VML 0)

#include_directories(${MKL_HOME}/include)
set(LAPACK_LIBRARY -L/home/00504/tg457645/lib -llapack -lblas -lgfortran)
set(BLAS_LIBRARY   -L/home/00504/tg457645/lib -lblas)

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

