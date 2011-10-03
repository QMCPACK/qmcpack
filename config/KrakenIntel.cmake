# tool chain for Kraken, Cray XT5  using PrgEnv-gnu
# switch the prog. env to use GNU compilers
# module swap PrgEnv-pgi PrgEnv-gnu

SET(CMAKE_SYSTEM_PROCESSOR "XT5")
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/4.9/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/4.9/bin/CC)

set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
#set(GNU_FLAGS "-fopenmp -O3 -ftemplate-depth-60 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
#set(XT_FLAGS "-march=amdfam10 -msse3 -D_CRAYMPI")
#set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
#set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS}")

set(INTEL_OPTS "-g -unroll -O3 -ip -openmp -opt-prefetch -ftz -xSSE3")
set(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${GNU_OPTS} ${INTEL_OPTS} -restrict -Wno-deprecated ")
set(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${INTEL_OPTS} -std=c99 -restrict -Wno-deprecated")

SET(CMAKE_Fortran_FLAGS "${INTEL_OPTS}")
SET(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
SET(USE_PREFETCH 1)
SET(PREFETCH_AHEAD 12)

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)

set(ACML_HOME /opt/acml/4.4.0/gfortran64/)

set(CMAKE_FIND_ROOT_PATH
      /opt/fftw/3.2.2.1
      /sw/xt/hdf5/1.8.5/cnl2.2_gnu4.4.3
      /sw/xt/szip/2.1/sles10.1_gnu4.4.3
      /nics/a/proj/qmc/kraken/intel/einspline
      /nics/a/proj/qmc/boost_1_38_0
     )

set(MKLROOT /opt/intel/Compiler/11.1/038/mkl)
include_directories(${MKLROOT}/include)
set(LAPACK_LIBRARY -Wl,--start-group ${MKLROOT}/lib/em64t/libmkl_intel_lp64.a ${MKLROOT}/lib/em64t/libmkl_sequential.a ${MKLROOT}/lib/em64t/libmkl_core.a -Wl,--end-group -lpthread)
