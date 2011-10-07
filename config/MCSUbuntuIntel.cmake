#--------------------------------------------------------------------------
# toolchain for Linux Clusters
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER  mpicc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(INTEL_OPTS "-g  -restrict -unroll  -O3 -openmp -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS} -std=c99")
set(MKL_HOME "/soft/com/packages/intel/12.0/174/mkl" CACHE STRING "MKL HOME")

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH
    /sandbox/libxml2-2.7.8
    /sandbox/boost_1_45_0
    /sandbox/einspline-0.9.2
    /sandbox/hdf5-1.6.6
    /sandbox/fftw3
    )

#--------------------------------------------------------------------------
# below is common for INTEL compilers and MKL library
#--------------------------------------------------------------------------
set(QMC_BUILD_STATIC 1)
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
set(LAPACK_LIBRARY -L${MKL_HOME}/lib/intel64 -lmkl_intel_lp64)
set(BLAS_LIBRARY -lmkl_sequential -lmkl_core -lpthread)
# set(LAPACK_LIBRARY /usr/lib/liblapack-3.a)
# set(BLAS_LIBRARY   /usr/lib/libblas-3.a)
# set(FORTRAN_LIBRARIES /usr/lib/gcc/x86_64-linux-gnu/4.4/libgfortran.a)
set(BUILD_SANDBOX 1)

INCLUDE(Platform/UnixPaths)

# FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
#   SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
#   SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
#   SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
#   SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
# ENDFOREACH(type)
# SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
# SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
# SET(CMAKE_C_LINK_SHARED_LIBRARY)
# SET(CMAKE_C_LINK_MODULE_LIBRARY)

