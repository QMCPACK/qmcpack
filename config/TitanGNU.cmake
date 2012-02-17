SET(CMAKE_SYSTEM_PROCESSOR "XK6")
#2012-02-16

set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/5.04/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/5.04/bin/CC)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS " -fomit-frame-pointer -malign-double  -fopenmp -O3 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS " -msse -msse2 -msse3 -D_CRAYMPI")
#4.6 is not working with the rest of the modules
#set(XT_FLAGS "-march=bdver1 -msse3 -D_CRAYMPI")
#set(XT_FLAGS "-march=amdfam10 -msse3 -D_CRAYMPI")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -std=c99")# -DHAVE_SSE -DHAVE_SSE2 ")

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
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
  /opt/cray/hdf5/1.8.6/gnu/45
  /opt/fftw/3.3.0.0/interlagos
  /ccs/home/jnkim/xk6/gnu45/libxml2
  /sw/xk6/boost/1.44.0/cle4.0_gnu4.5.3
)

# bypass einspline search and build it inside QMCPACK
set(EINSPLINE_HOME /lustre/widow3/scratch/jnkim/einspline)
set(HAVE_EINSPLINE 1)
set(HAVE_EINSPLINE_EXT 0)

#link_libraries(/opt/xt-libsci/11.0.04.4/gnu/45/interlagos/lib/libsci_gnu.a)
