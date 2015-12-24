#--------------------------------------------------------------------------
# toolchain for Linux Clusters
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# setting compilers, compiler options and MKL_HOME
#--------------------------------------------------------------------------
set(CMAKE_CXX_COMPILER /opt/apps/gnu44/openmpi-1.6.2/bin/mpiCC)
set(CMAKE_C_COMPILER  /opt/apps/gnu44/openmpi-1.6.2/bin/mpicc)
#set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER -DUSE_REAL_STRUCT_FACTOR")
#set(INTEL_OPTS "-g  -restrict -unroll  -O3 -ip  -xSSE4.2 -openmp -Wno-deprecated")
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER=1 -DUSE_REAL_STRUCT_FACTOR")
set(GNU_FLAGS "-malign-double -fomit-frame-pointer -ffast-math -fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-march=amdfam10 -msse3")
#set(XT_FLAGS "-march=bdver1 -D_CRAYMPI -DHAVE_FMA4=1 -DHAVE_AMDLIBM=1")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -std=c99")

#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH
 /opt/apps/gnu44/fftw-3.3.2
 /opt/apps/gnu44/hdf5-1.8.10
 /opt/apps/gnu44/openmpi-1.6.2
 /usr/lib64
)

#--------------------------------------------------------------------------
# below is common for INTEL compilers and MKL library
#--------------------------------------------------------------------------
set(ENABLE_OPENMP 1)
SET(QMC_BUILD_STATIC 1)
set(HAVE_MPI 1)
set(HAVE_SSE 1)
set(HAVE_SSE2 1)
set(HAVE_SSE3 1)
set(HAVE_SSSE3 1)
set(HAVE_SSE41 0)
set(USE_PREFETCH 1)
set(PREFETCH_AHEAD 10)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_SHARED_LINKER_FLAGS "")

include_directories(/usr/include/libxml2)
link_libraries(/usr/lib64/libxml2.so)
link_libraries(/usr/lib64/libgfortran.so.3.0.0)
link_libraries(/usr/lib64/liblapack.so)
