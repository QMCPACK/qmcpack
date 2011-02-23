SET(CMAKE_SYSTEM_PROCESSOR "XT5")


set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/4.7/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/4.7/bin/CC)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
#set(GNU_FLAGS "-Wl,-z,muldefs -fopenmp -O3 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(GNU_FLAGS "-fopenmp -O3 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-march=amdfam10 -msse3 -D_CRAYMPI")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS}")

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)


set(CMAKE_FIND_ROOT_PATH
  /global/homes/j/jkim/share/hopper/einspline/gnu_4.5.2
  /opt/cray/hdf5/1.8.5.0/hdf5-gnu
  /opt/fftw/3.2.2.1
  /global/homes/j/jkim/share/boost_1_42_0
  /global/homes/j/jkim/share/hopper/libxml2/cnos_gnu_4.5.2
)

SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
SET(USE_PREFETCH 1)
SET(PREFETCH_AHEAD 12)

set(ACML_HOME /opt/acml/4.4.0/gnu64)
SET(ACML_LIBRARIES ${ACML_HOME}/lib/libacml.a ${ACML_HOME}/lib/libacml_mv.a)
link_libraries(${ACML_LIBRARIES})
