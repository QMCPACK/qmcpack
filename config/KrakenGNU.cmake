# tool chain for Kraken, Cray XT5  using PrgEnv-gnu
# - switch the prog. env to use GNU compilers PrgEnv-gnu
# - remove acml
# Here is the module used by jnkim on 2012-09-06
#  1) modules/3.2.6.6
#  2) portals/2.2.0-1.0301.26633.6.9.ss
#  3) nodestat/2.2-1.0301.25918.4.1.ss
#  4) sdb/1.0-1.0301.25929.4.88.ss
#  5) MySQL/5.0.64-1.0301.2899.20.1.ss
#  6) lustre-cray_ss_s/1.8.4_2.6.27.48_0.12.1_1.0301.5943.18.1-1.0301.27524.1.24
#  7) Base-opts/1.0.2-1.0301.25878.4.1.ss
#  8) xtpe-network-seastar
#  9) torque/2.5.11
# 10) moab/6.1.6
# 11) xtpe-istanbul
# 12) gold
# 13) hsi
# 14) tgusage/3.0-r3
# 15) altd/1.0
# 16) globus/5.0.4
# 17) DefApps
# 18) gcc/4.6.2
# 19) xt-libsci/11.0.04
# 20) xt-mpich2/5.3.5
# 21) pmi/2.1.4-1.0000.8596.15.1.ss
# 22) xt-asyncpe/5.11
# 23) atp/1.4.1
# 24) PrgEnv-gnu/3.1.72
# 25) subversion/1.6.9
# 26) numpy/1.6.1
# 27) cmake/2.8.2
# 28) hdf5/1.8.6
# 29) fftw/3.3.0.0


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

set(ACML_HOME /opt/acml/4.4.0/gfortran64/)

set(CMAKE_FIND_ROOT_PATH
  /opt/cray/hdf5/1.8.6/gnu/46
  /opt/fftw/3.3.0.0/x86_64/
  /nics/a/proj/qmc/boost_1_38_0
  /nics/a/proj/qmc
  )

#uncomment this if external build to be used
#set(HAVE_EINSPLINE_EXT 1)
