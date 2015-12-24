#--------------------------------------------------------------------------
# Created by Raymond Clay on 10/2/2012 for Taub Campus Cluster
# 
#Currently Loaded Modulefiles:
#  1) torque/3.0.6        3) env/taub            5) intel/11.1          7) svn/1.6             9) boost/1.51.0
#  2) moab/6.1.8          4) vim/7.3             6) openmpi/1.4-intel   8) cmake/2.8

#--------------------------------------------------------------------------
SET(CMAKE_SYSTEM_PROCESSOR "x86_64")

#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER mpicc)
SET(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
SET(INTEL_OPTS "-g -std=c99 -restrict -unroll  -O3 -ip -xT -openmp -Wno-deprecated")
SET(MKL_HOME "/usr/local/intel-11.1/mkl")
SET(MKL_EXT "em64t")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS}")

SET(CMAKE_Fortran_FLAGS "${INTEL_OPTS}")
SET(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH
  /home/rcclay2
  /usr/local/boost/boost-1.51.0
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

include_directories(/usr/local/intel-11.1/mkl/include)
link_libraries(-L/usr/local/intel-11.1/mkl/lib/em64t -mkl=sequential)

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

