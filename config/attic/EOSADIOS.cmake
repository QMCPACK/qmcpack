#######################################
# module load                         #
######################################
SET(CMAKE_SYSTEM_PROCESSOR "XC30")

set(HAVE_ADIOS 1)

set(CMAKE_CXX_COMPILER  /opt/cray/craype/default/bin/CC)
set(CMAKE_C_COMPILER  /opt/cray/craype/default/bin/cc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER -DUSE_REAL_STRUCT_FACTOR")
set(INTEL_OPTS " -static -g  -restrict -unroll  -O3 -ip  -xAVX -openmp -Wno-deprecated")
set(XT_FLAGS "-D_CRAYMPI")# -DHAVE_FMA4=1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XT_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${XT_FLAGS} ${INTEL_OPTS} -std=c99")

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
set(HAVE_MKL 1)
set(HAVE_MKL_VML 1)

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
  /opt/fftw/3.3.0.1/x86_64
  /opt/cray/hdf5/1.8.11/INTEL/130
  /opt/fftw/3.3.0.4/sandybridge
  /ccs/home/jnkim/eos/libxml2
  /sw/xc30/boost/1.54.0/cle4_intel13.1.3.192/
)
include_directories($ENV{MKLROOT}/include)
link_libraries(-L$ENV{MKLROOT}/lib/intel64 -mkl=sequential)

#add ADIOS stuff provided by adios_config 
#DIR=/sw/xc30/adios/1.5.0/sles11.2_intel13.1
#CFLAGS=-I/sw/xc30/adios/1.5.0/sles11.2_intel13.1/include -I/sw/xc30/mxml/2.6/sles11.2_intel13.1/include
#LDFLAGS=-L/sw/xc30/adios/1.5.0/sles11.2_intel13.1/lib -ladios -L/sw/xc30/mxml/2.6/sles11.2_intel13.1/lib -lmxml -lm -lmxml
include_directories(
  /sw/xc30/adios/1.5.0/sles11.2_intel13.1/include
  /sw/xc30/mxml/2.6/sles11.2_intel13.1/include
)
link_libraries(-L/sw/xc30/adios/1.5.0/sles11.2_intel13.1/lib -ladios -L/sw/xc30/mxml/2.6/sles11.2_intel13.1/lib -lmxml -lm -lmxml)
