# tool chain for Kraken, Cray XT5  using PrgEnv-gnu
#Currently Loaded Modulefiles:
#  1) modules/3.2.6.6
#  2) portals/2.2.0-1.0301.26633.6.9.ss
#  3) nodestat/2.2-1.0301.25918.4.1.ss
#  4) sdb/1.0-1.0301.25929.4.88.ss
#  5) MySQL/5.0.64-1.0301.2899.20.1.ss
#  6) lustre-cray_ss_s/1.8.4_2.6.27.48_0.12.1_1.0301.5943.18.1-1.0301.27524.1.24
#  7) Base-opts/1.0.2-1.0301.25878.4.1.ss
#  8) xtpe-network-seastar
#  9) PrgEnv-gnu/3.1.72
# 10) atp/1.4.1
# 11) xt-asyncpe/5.11
# 12) pmi/2.1.4-1.0000.8596.15.1.ss
# 13) xt-mpich2/5.3.5
# 14) xt-libsci/11.0.04
# 15) gcc/4.6.2
# 16) torque/2.5.11
# 17) moab/6.1.6
# 18) xtpe-istanbul
# 19) gold
# 20) hsi
# 21) tgusage/3.0-r3
# 22) altd/1.0
# 23) globus/5.0.4
# 24) xdusage/1.0-r7
# 25) DefApps
# 26) hdf5/1.8.6
# 27) fftw/3.3.0.0
# 28) cmake/2.8.7
# 29) boost/1.44.0
# 30) subversion/1.6.9


SET(CMAKE_SYSTEM_PROCESSOR "XT5")
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/default/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/default/bin/CC)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS "-fopenmp -O3 -ftemplate-depth-60 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-march=amdfam10 -msse2 -D_CRAYMPI")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -std=c99")

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

#set(ACML_HOME /opt/acml/4.4.0/gfortran64/)

set(CMAKE_FIND_ROOT_PATH
  /opt/cray/hdf5/1.8.6/gnu/46
  /opt/fftw/3.3.0.0/x86_64/
  /sw/xt-cle3.1/boost/1.44.0/cnl3.1_gnu4.6.1/include
  )

set(BOOST_HOME /sw/xt-cle3.1/boost/1.44.0/cnl3.1_gnu4.6.1)

#uncomment this if external build to be used
#set(HAVE_EINSPLINE_EXT 1)
