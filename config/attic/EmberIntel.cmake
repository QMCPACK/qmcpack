#--------------------------------------------------------------------------
# tool chain for ember Wed Sep 24 2010
# Currently Loaded Modulefiles:
#  1) sgi-mpt-2.01         3) mkl-11.1             5) mssftp-client-2.9    7) saveafterjob-2.1.0   9) globus-4.0.8-r2     11) tgusage-3.0
#  2) intel-11.1.073       4) cue-mkl              6) uberftp-client-2.6   8) gx-map-0.5.3.3-r1   10) java-sun-1.6        12) tg-policy-0.2-r1
#--------------------------------------------------------------------------
SET(CMAKE_SYSTEM_PROCESSOR "EMBER")

#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER  icpc)
set(CMAKE_C_COMPILER  icc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")

#set(INTEL_OPTS "-g -unroll -ansi-alias -O3 -ip -openmp -opt-prefetch -ftz -xSSE4.2")
set(INTEL_OPTS "-g -unroll -O3 -ip -openmp -opt-prefetch -ftz -xSSE4.2")
set(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${GNU_OPTS} ${INTEL_OPTS} -restrict -Wno-deprecated")
set(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${INTEL_OPTS} -std=c99 -restrict -Wno-deprecated")

SET(CMAKE_Fortran_FLAGS "${INTEL_OPTS}")
SET(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH
  /u/ac/jnkim/share/gnu/hdf5
  /u/ac/jnkim/share/intel/einspline_sse
  /u/ac/jnkim/share/boost
  /usr/apps/math/fftw/intel/3.2.1
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

set(MKL_HOME "/usr/local/intel/Compiler/11.1/072/mkl")
include_directories(${MKL_HOME}/include)
link_libraries(-L/usr/local/intel/Compiler/11.1/072/mkl/lib/em64t -mkl=sequential -lmpi)

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

