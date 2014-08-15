SET(CMAKE_SYSTEM_PROCESSOR "XK7")
#2011-12-06

set(HAVE_ADIOS 1)
set(CMAKE_C_COMPILER  cc)
set(CMAKE_CXX_COMPILER  CC)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER=1 -DUSE_REAL_STRUCT_FACTOR") 
#set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DUSE_REAL_STRUCT_FACTOR -DDISABLE_TIMER=1 -DHAVE_FMA4=1 -DHAVE_AMDLIBM=1")
set(GNU_FLAGS "-malign-double -fomit-frame-pointer -ffast-math -fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated")
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
  /opt/cray/hdf5/1.8.8/gnu/47
  /opt/fftw/3.3.0.3/x86_64/lib
  /sw/xk6/boost/1.44.0/cle4.0_gnu4.5.3
  /ccs/home/jnkim/titan/gnu47/libxml2
)

#AMD math lib
include_directories(/ccs/home/jnkim/lib/amdlibm/include)
link_libraries(/ccs/home/jnkim/lib/amdlibm/lib/static/libamdlibm.a)

include_directories(
  /ccs/proj/e2e/pnorbert/ADIOS/chester.gnu/include
  #/opt/cray/pmi/default/include 
  #/opt/cray/gni-headers/default/include
  #/sw/xk6/adios/1.6.0/cle4.0_gnu4.7.2/include 
  #/sw/xk6/mxml/2.6/cle4.0_gnu4.5.3/include 
  #/sw/xk6/dataspaces/1.3.0/cle4.0_gnu4.7.2/include 
  )
link_libraries(
-L/ccs/proj/e2e/pnorbert/ADIOS/chester.gnu/lib -ladios -L/sw/xk6/mxml/2.6/cle4.0_gnu4.5.3/lib -lmxml -L/sw/xk6/dataspaces/1.4.0/cle4.2_gnu4.8.2/lib -L/sw/xk6/dataspaces/1.4.0/cle4.2_gnu4.8.2/lib -L/ccs/proj/e2e/chaos/titan/gnu/lib -L/ccs/proj/e2e/chaos/titan/gnu/lib -L/ccs/proj/e2e/chaos/titan/gnu/lib -L/ccs/proj/e2e/chaos/titan/gnu/lib -L/ccs/proj/e2e/chaos/titan/gnu/lib -L/opt/cray/pmi/default/lib64 -L/opt/cray/ugni/default/lib -lm -lmxml -ldspaces -ldscommon -ldart -ldspaces -ldscommon -ldart -levpath -lffs -latl -ldill -lcercs_env -lpmi -lugni
  #-L/sw/xk6/adios/1.6.0/cle4.0_gnu4.7.2/lib -ladios -L/sw/xk6/mxml/2.6/cle4.0_gnu4.5.3/lib -lmxml -L/sw/xk6/dataspaces/1.3.0/cle4.0_gnu4.7.2/lib -L/sw/xk6/dataspaces/1.3.0/cle4.0_gnu4.7.2/lib -L/opt/cray/pmi/default/lib64 -L/opt/cray/ugni/default/lib64 -lm -lmxml -ldspaces -ldscommon -ldart -ldspaces -ldscommon -ldart -lpmi -lugni
  #-L/ccs/home/zgu/adioshub/adios_titan/lib -ladios -lm
  #-L/sw/xk6/mxml/2.6/cle4.0_gnu4.5.3/lib -lmxml
  -lz
  )
