#--------------------------------------------------------------------------
# toolchain for keeneland cuda-enabled 
# cuda is broken with intel compilers : suspect cmake's cuda module
# important to use PE-gnu environment to use openmpi/1.4.3-gnu
# module swap PE-intel PE-gnu
# Currently Loaded Modulefiles:
#1) modules             3) moab/5.4.0.s16449   5) hsi                 7) openmpi/1.5.1-gnu   9) automake/1.9.6     11) subversion/1.6.15
#2) torque/2.4.8        4) gold                6) PE-gnu              8) cuda/3.2           10) autoconf/2.68      12) cmake/2.8.0
#--------------------------------------------------------------------------
# setting compilers, compiler options 
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER /sw/keeneland/openmpi/1.5.1/centos5.5_gnu4.1.2/bin/mpicxx)
set(CMAKE_C_COMPILER gcc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")

#set(INTEL_OPTS "-g -unroll -ansi-alias -O3 -ip -openmp -opt-prefetch -ftz -xSSE4.2")
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS "-fopenmp -O3 -ftemplate-depth-60 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-march=nocona -msse3 ")

set(CMAKE_CXX_FLAGS "-g ${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "-g ${XT_FLAGS} ${GNU_FLAGS}")

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH
  /nics/a/proj/qmc/keeneland/gnu4.1/einspline
  /sw/keeneland/fftw/3.2.1/centos5.4_gnu4.1.2_fPIC
  /sw/keeneland/hdf5/1.8.6/centos5.5_gnu4.1.2
  /sw/keeneland/szip/2.1/centos5.5_gnu4.1.2
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

link_libraries(-llapack -lblas)
INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)
#
