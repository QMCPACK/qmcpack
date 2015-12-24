#--------------------------------------------------------------------------
# tool chain for ember Wed Sep 24 2010
# Currently Loaded Modulefiles:
#Currently Loaded Modulefiles:
#1) modules               4) gold                  7) intel/2011.1.107     10) cuda/3.2             13) subversion/1.6.15
#2) torque/2.4.8          5) hsi                   8) openmpi/1.5.1-intel  11) automake/1.9.6       14) cmake/2.8.0
#3) moab/5.4.0.s16449     6) mkl/2011.1.107        9) PE-intel             12) autoconf/2.68
#--------------------------------------------------------------------------
SET(CMAKE_SYSTEM_PROCESSOR "EMBER")

#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER /usr/local/mpi/mvapich2/mvapich2-1.7rc1-intel-12.0.4/bin/mpicxx)
set(CMAKE_C_COMPILER  icc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")

#set(INTEL_OPTS "-g -unroll -ansi-alias -O3 -ip -openmp -opt-prefetch -ftz -xSSE3")
set(INTEL_OPTS "-g -unroll -O3 -ip -openmp -opt-prefetch -ftz -xSSE2")
set(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${GNU_OPTS} ${INTEL_OPTS} -restrict -Wno-deprecated")
set(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${INTEL_OPTS} -std=c99 -restrict -Wno-deprecated")

SET(CMAKE_Fortran_FLAGS "${INTEL_OPTS}")
SET(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH
		/usr/local/hdf/phdf5/v187
		/usr/apps/math/fftw/intel/3.3.0/mvapich2-1.7rc1
		/usr/local/hdf/szip/v2.1/static/encoder
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

include_directories(/usr/local/intel/mkl/include)
link_libraries(-L/usr/local/intel/mkl/lib/intel64 -mkl=sequential)

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

