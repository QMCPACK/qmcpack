#######################################
# toolchain files on Macbook Pro
# 2015-10-10
# OS X 10.10.5 
# Intel XE 15.0.x and higher
# Jeongnim Kim, Intel
#######################################
# MPICH compiled with Intel compilers
set(CMAKE_CXX_COMPILER  /usr/local/bin/mpicxx)
set(CMAKE_C_COMPILER  icc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER -DUSE_REAL_STRUCT_FACTOR")
set(INTEL_OPTS " -g  -restrict -unroll  -O3 -ip  -xHOST -openmp -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS} -std=c99")

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

set(FFTW_INCLUDE_DIR $ENV{MKLROOT}/include/fftw)
set(FFTW_LIBRARIES "")

include_directories($ENV{MKLROOT}/include $ENV{MKLROOT}/include/fftw)
link_libraries(-mkl)

