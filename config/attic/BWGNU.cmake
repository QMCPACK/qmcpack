SET(CMAKE_SYSTEM_PROCESSOR "XK6")
#2011-12-06

set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/default/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/default/bin/CC)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
#set(GNU_FLAGS "-fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated -ffast-math -funroll-loops -fomit-frame-pointer -malign-double")
set(GNU_FLAGS "-fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated -ffast-math -funroll-loops -fomit-frame-pointer -malign-double -DDISABLE_TIMER=1 -DHAVE_FMA4=1 -DHAVE_AMDLIBM=1")
#set(GNU_FLAGS "-fopenmp -g -O0 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated -ffast-math -funroll-loops -fomit-frame-pointer -malign-double")
#set(XT_FLAGS "-msse3 -D_CRAYMPI") 
#set(XT_FLAGS "-march=amdfam10 -msse3 -D_CRAYMPI")
#interlogs bdver1 but without it better
set(XT_FLAGS "-march=bdver1  -D_CRAYMPI") 
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -std=c99")

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
SET(HAVE_SSE41 1)
SET(USE_PREFETCH 1)
SET(PREFETCH_AHEAD 12)
#SET(HAVE_AMDLIBM 1)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_SHARED_LINKER_FLAGS "")

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-static")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-static")
ENDFOREACH(type)

set(CMAKE_FIND_ROOT_PATH
/opt/cray/hdf5/default/gnu/47
/opt/fftw/default/interlagos
/u/staff/rmokos/libs/boost
/u/staff/rmokos/libs/libxml2
)

#set(HAVE_LIBBOOST 1)
#include_directories(/u/staff/rmokos/libs/boost)

#set(EINSPLINE_HOME /u/sciteam/jnkim/svnwork/einspline)
#set(HAVE_EINSPLINE 1)
#set(EINSPLINE_SSE_BUG 1)
#set(HAVE_EINSPLINE_EXT 0)
#link_libraries(/opt/acml/5.0.0/gfortran64_int64/lib/libacml.a)

# for Torsten's libPGT
#include_directories(/u/staff/rmokos/libs/libpgt/gnu/include/libpgt-0.2)
#link_libraries(/u/staff/rmokos/libs/libpgt/gnu/lib/libpgt-0.2.a)

# for AMD's libm
include_directories(/u/staff/rmokos/libs/amdlibm/include)
link_libraries(/u/staff/rmokos/libs/amdlibm/lib/static/libamdlibm.a)

