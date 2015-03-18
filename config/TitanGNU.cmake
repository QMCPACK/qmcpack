SET(CMAKE_SYSTEM_PROCESSOR "XK7")
# Last updated: Mar 16, 2015
# NEED THESES + defaults
#  module swap PrgEnv-pgi PrgEnv-gnu
#  module load cray-hdf5
#  module load fftw
#  module load boost
#  module load subversion
#  module load cmake

set(CMAKE_C_COMPILER  cc)
set(CMAKE_CXX_COMPILER  CC)

set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER=1 -DUSE_REAL_STRUCT_FACTOR") 
#set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DUSE_REAL_STRUCT_FACTOR -DDISABLE_TIMER=1 -DHAVE_FMA4=1 -DHAVE_AMDLIBM=1")
set(GNU_FLAGS "-malign-double -fomit-frame-pointer -ffast-math -fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-march=bdver1 -D_CRAYMPI -DHAVE_FMA4=1 -DHAVE_AMDLIBM=1")
#set(XT_FLAGS "-msse3 -D_CRAYMPI")
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
  $ENV{HDF5_DIR}
  $ENV{FFTW_DIR}/..
  $ENV{BOOST_DIR}
  /sw/xk7/libxml2
)

#AMD math lib
include_directories(/sw/xk7/amdlibm/include)
link_libraries(/sw/xk7/amdlibm/lib/static/libamdlibm.a)

