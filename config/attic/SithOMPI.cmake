
######################################
#    modules load                    #
######################################
#source /etc/profile.d/modules.sh
#module unload pgi PE-pgi
#module load PE-gnu/4.7.1
#module load ompi/1.6.5
#module load netcdf/4.1.3
#module load mxml/2.7
#module load adios/1.5.0
#module load dataspaces
#module load acml/5.3.0
#module load szip/2.1
#module load hdf5-parallel/1.8.11
#module load cmake/2.8.10.2
#module load subversion

#export HDF5_HOME=/sw/redhat6/hdf5-parallel/1.8.11/rhel6.4_gnu4.7.1

SET(CMAKE_SYSTEM_PROCESSOR "XK7")

set(CMAKE_C_COMPILER  gcc)
set(CMAKE_CXX_COMPILER  /sw/sith/ompi/1.6.5/rhel6_gnu4.7.1/bin/mpicxx)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER=1 -DUSE_REAL_STRUCT_FACTOR") 
set(GNU_FLAGS "-malign-double -fomit-frame-pointer -ffast-math -fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
set(XT_FLAGS "-march=amdfam10 -msse3")
#set(XT_FLAGS "-march=bdver1 -D_CRAYMPI -DHAVE_FMA4=1 -DHAVE_AMDLIBM=1")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -std=c99")

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

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_SHARED_LINKER_FLAGS "")

#FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
#  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
#  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-static")
#  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
#  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-static")
#ENDFOREACH(type)

set(CMAKE_FIND_ROOT_PATH
  /sw/redhat6/hdf5/1.8.11/rhel6.4_gnu4.7.1
  /sw/redhat6/szip/2.1/rhel6_gnu4.7.1
  /ccs/home/zgu/Software_sith
)

link_libraries(/sw/redhat6/acml/5.3.0/gfortran64/lib/libacml.a -lz)
link_libraries(/ccs/compilers/gcc/rhel6-x86_64/4.7.1/lib64/libgfortran.so)
