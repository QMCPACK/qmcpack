#--------------------------------------------------------------------------
# tool chain for abe Wed Aug 12 2009
#--------------------------------------------------------------------------
SET(CMAKE_SYSTEM_PROCESSOR "ES")

#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
#set(OPENMPI_HOME /home/jnkim/share/intel/openmpi)
set(CMAKE_CXX_COMPILER  /opt/intel/impi/4.0.0/intel64/bin/mpiicpc)
set(CMAKE_C_COMPILER  icc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
#set(INTEL_OPTS "-g -unroll -ansi-alias -O3 -ip -openmp -opt-prefetch -ftz -xSSE4.2")
set(INTEL_OPTS "-g -unroll -O3 -ip -openmp -opt-prefetch -ftz -xSSE4.2")
#SET(CMAKE_CXX_FLAGS "-g -O3 -msse4.2 -ftemplate-depth-60 -Drestrict=__restrict__ -funroll-all-loops   -finline-limit=1000 -Wno-deprecated ")
#SET(CMAKE_C_FLAGS "-O3 -msse4.2 -Drestrict=__restrict__ -funroll-all-loops   -finline-limit=1000 -std=gnu99 -fomit-frame-pointer ")
set(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${GNU_OPTS} ${INTEL_OPTS} -restrict -Wno-deprecated")
set(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${INTEL_OPTS} -std=c99 -restrict -Wno-deprecated")
SET(CMAKE_Fortran_FLAGS "${INTEL_OPTS}")
SET(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS} -std=c99")

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH
    /home/jnkim/share/boost
    /home/jnkim/share/intel/einspline
    /home/jnkim/share/intel/fftw-3.2s
    /home/jnkim/share
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

set(MKL_HOME "/usr/local/intel/11.1.056/mkl")
include_directories(${MKL_HOME}/include)
#include_directories(/usr/lib64/mpi/gcc/openmpi/include)
#link_libraries(-L/usr/local/intel/11.1.056/mkl/lib/em64t  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core) #-lm -lpthread)
link_libraries(-L/usr/local/intel/11.1.056/mkl/lib/em64t -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread) 
#-L${OPENMPI_HOME}/lib -lmpi -lopen-rte -lopen-pal -ldl)


INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

