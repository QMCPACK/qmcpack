SET(CMAKE_SYSTEM_PROCESSOR "XK6")
#2011-12-06
#Currently Loaded Modulefiles:
#  1) modules/3.2.6.6
#  2) DefApps
#  3) torque/2.4.1b1-snap.200905191614
#  4) moab/5.3.6
#  5) nodestat/2.2-1.0400.29866.4.3.gem
#  6) sdb/1.0-1.0400.30000.6.18.gem
#  7) MySQL/5.0.64-1.0000.4667.20.1
#  8) lustre-cray_gem_s/1.8.4_2.6.32.45_0.3.2_1.0400.6221.1.1-1.0400.30303.0.0novmap1
#  9) udreg/2.3.1-1.0400.3911.5.6.gem
# 10) ugni/2.3-1.0400.3912.4.29.gem
# 11) gni-headers/2.1-1.0400.3906.5.1.gem
# 12) dmapp/3.2.1-1.0400.3965.10.12.gem
# 13) xpmem/0.1-2.0400.29883.4.6.gem
# 14) hss-llm/6.0.0
# 15) Base-opts/1.0.2-1.0400.29823.8.1.gem
# 16) xtpe-network-gemini
# 17) PrgEnv-gnu/4.0.30
# 18) xt-mpich2/5.4.0
# 19) atp/1.4.0
# 20) xt-asyncpe/5.04
# 21) pmi/3.0.0-1.0000.8661.28.2807.gem
# 22) xt-libsci/11.0.04.4
# 23) gcc/4.5.3
# 24) xtpe-interlagos
# 25) subversion/1.6.17
# 26) hdf5/1.8.6
# 27) fftw/3.3.0.0
# 28) cmake/2.8.6
# 29) boost/1.44.0


set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/5.04/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/5.04/bin/CC)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS "-fopenmp -O3 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
#set(XT_FLAGS "-march=amdfam10 -msse3 -D_CRAYMPI")
set(XT_FLAGS "-msse3 -D_CRAYMPI")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS}")

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
  /ccs/home/jnkim/xk6/gnu45
  /sw/xk6/boost/1.44.0/cle4.0_gnu4.5.3
)

#link_libraries(/opt/xt-libsci/11.0.04.4/gnu/45/interlagos/lib/libsci_gnu.a)
