#--------------------------------------------------------------------------
# toolchain for Linux Clusters
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER  icc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(INTEL_OPTS "-g  -restrict -unroll  -O3 -ip -xT -openmp -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS} -std=c99")
set(MKL_HOME "/opt/intel/mkl/10.0.3.020" CACHE STRING "MKL HOME")

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH
    $ENV{EINSPLINE_HOME}
    $ENV{BOOST_HOME}
    $ENV{HDF5_HOME}
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
set(HAVE_MKL 1)
set(HAVE_MKL_VML 1)

include_directories(${MKL_HOME}/include)
set(LAPACK_LIBRARY -L${MKL_HOME}/lib/em64t -lmkl_lapack)
set(BLAS_LIBRARY -L${MKL_HOME}/lib/em64t -lmkl -lguide)

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

