#--------------------------------------------------------------------------
# tool chain for abe Wed Aug 12 2009
#--------------------------------------------------------------------------
SET(CMAKE_SYSTEM_PROCESSOR "ES")

#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER /usr/local/mvapich2-1.2-intel-ofed-1.2.5.5/bin/mpicxx)
set(CMAKE_C_COMPILER  /usr/local/intel/10.1.017/bin/icc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -Drestrict=__restrict__")
set(INTEL_OPTS "-g  -restrict -unroll  -O3 -ip -xT -openmp -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS} -std=c99")
set(MKL_HOME "/usr/local/intel/mkl/10.1.2.024" CACHE STRING "MKL HOME")
set(NVCC_CXX_FLAGS "-arch;sm_13;-Drestrict=__restrict")

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(EINSPLINE_HOME "/u/ac/esler/svnwork/einspline" CACHE STRING "Einspline source directory")
set(CMAKE_FIND_ROOT_PATH
    /u/ac/jnkim/share/boost
    /usr/apps/hdf/hdf5/v167
    /usr/apps/hdf/szip/v2.1/static/encoder
    /usr/local/libxml2-2.6.29
    /usr/apps/math/fftw/fftw-3.1.2/intel10
    )

#--------------------------------------------------------------------------
# below is common for INTEL compilers and MKL library
#--------------------------------------------------------------------------
set(ENABLE_OPENMP 1)
set(QMC_CUDA 1)
set(HAVE_MPI 1)
set(HAVE_SSE 1)
set(HAVE_SSE2 1)
set(HAVE_SSE3 1)
set(HAVE_SSSE3 1)
set(USE_PREFETCH 1)
set(PREFETCH_AHEAD 10)
set(HAVE_MKL 1)
set(HAVE_MKL_VML 1)

ADD_DEFINITIONS(-Drestrict=__restrict__)

include_directories(${MKL_HOME}/include)
set(LAPACK_LIBRARY -L${MKL_HOME}/lib/em64t -lmkl_lapack)
set(BLAS_LIBRARY -L${MKL_HOME}/lib/em64t -lmkl -lguide)

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

