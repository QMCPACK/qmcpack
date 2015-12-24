#--------------------------------------------------------------------------
# toolchain for Intel Phi on beacon using native mode
# only build einspline 
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_C_COMPILER  icc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
#set(INTEL_OPTS "-g  -restrict -unroll  -O3 -ip -xT -openmp -Wno-deprecated")
set(INTEL_OPTS "-restrict -unroll  -O3 -ip -mmic -openmp -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS} ${GNU_OPTS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS} -std=c99")

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
#set(CMAKE_FIND_ROOT_PATH
#    $ENV{EINSPLINE_HOME}
#    $ENV{BOOST_HOME}
#    $ENV{HDF5_HOME}
#)

#--------------------------------------------------------------------------
# below is common for INTEL compilers and MKL library
#--------------------------------------------------------------------------
set(QMC_PHI 1)
set(ENABLE_OPENMP 1)
set(HAVE_MPI 0)
set(HAVE_SSE 0)
set(HAVE_SSE2 0)
set(HAVE_SSE3 0)
set(HAVE_SSSE3 0)
set(HAVE_SSE41 0)
set(USE_PREFETCH 0)
set(PREFETCH_AHEAD 10)
set(HAVE_MKL 1)
set(HAVE_MKL_VML 1)

#include_directories($ENV{MKLROOT}/include)
#link_libraries(-L$ENV{MKLROOT}/lib/intel64 -mkl=sequential)

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

