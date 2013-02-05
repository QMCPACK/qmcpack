SET(CMAKE_SYSTEM_PROCESSOR "XT5")
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

set(CMAKE_C_COMPILER  /opt/cray/xt-asyncpe/default/bin/cc)
set(CMAKE_CXX_COMPILER  /opt/cray/xt-asyncpe/default/bin/CC)
set(CMAKE_Fortran_COMPILER /opt/cray/xt-asyncpe/default/bin/ftn)

set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER=1  -DUSE_REAL_STRUCT_FACTOR")
#set(GNU_FLAGS "-Wl,-z,muldefs -fopenmp -O3 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(GNU_FLAGS "-malign-double -fomit-frame-pointer -ffast-math -fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-march=amdfam10 -msse3 -D_CRAYMPI")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -std=c99")
#set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS}")
#set(CMAKE_Fortran_FLAGS "-O3 -march=amdfam10 -fopenmp -funroll-all-loops")
#set(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})
##set(CMAKE_Fortran_FLAGS_DEBUG  "-march=amdfam10 -fopenmp  -msse3 -fno-f2c -O0 -g")
#set(CMAKE_Fortran_FLAGS_DEBUG  "-march=amdfam10 -fopenmp  -msse3 -O0 -g")

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)


set(CMAKE_FIND_ROOT_PATH
  /opt/cray/hdf5/1.8.8/gnu/47
  /opt/fftw/3.3.0.1/x86_64
  /usr/common/usg/boost/1.47.0
  /global/homes/j/jkim/share/hopper/libxml2/cnos_gnu_4.5.2
)
  #/global/homes/j/jkim/share/hopper/einspline/gnu_4.5.2

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
SET(HAVE_SSE41 0)
SET(USE_PREFETCH 1)
SET(PREFETCH_AHEAD 12)

#set(ACML_HOME /opt/acml/4.4.0/gnu64)
#SET(ACML_LIBRARIES ${ACML_HOME}/lib/libacml.a ${ACML_HOME}/lib/libacml_mv.a)
#link_libraries(${ACML_LIBRARIES})
