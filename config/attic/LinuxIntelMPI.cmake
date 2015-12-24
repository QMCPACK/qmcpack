#--------------------------------------------------------------------------
# toolchain for Linux Clusters
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER  icc)
#set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER -DUSE_REAL_STRUCT_FACTOR")
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER")
set(INTEL_OPTS "-g  -restrict -unroll  -O3 -ip  -xSSE4.2 -openmp -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS} -std=c99")

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH
  /usr/gapps/qmc/libs/INTEL/hdf5-1.8.8
  /usr/gapps/qmc/libs/INTEL/libxml2-2.7.4
  /usr/gapps/qmc/libs/INTEL/boost_1_40_0 
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
set(HAVE_SSE41 1)
set(USE_PREFETCH 1)
set(PREFETCH_AHEAD 10)
set(HAVE_MKL 1)
set(HAVE_MKL_VML 1)

include_directories($ENV{MKLROOT}/include)
link_libraries(-L$ENV{MKLROOT}/lib/intel64 -mkl=sequential)

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

