

######################################
#    modules load                    #
######################################
#source /etc/profile.d/modules.sh
#module unload pgi PE-pgi
#module load PE-gnu/4.7.1
#module load ompi/1.6.5
#module load netcdf/4.1.3
#module load mxml/2.7
#module load dataspaces
#module load acml/5.3.0
#module load szip/2.1
#module load hdf5-parallel/1.8.11
#module load cmake/2.8.10.2
#module load subversion

#export HDF5_HOME=/sw/redhat6/hdf5-parallel/1.8.11/rhel6.4_gnu4.7.1

SET(CMAKE_SYSTEM_PROCESSOR "XK7")

#set(HAVE_ADIOS 0)

#set(ADIOS_HOME /ccs/home/zgu/adioshub/adios_sith)
set(ADIOS_HOME /ccs/proj/e2e/pnorbert/ADIOS/sith.gnu)
set(ADIOS_INCLUDES ${ADIOS_HOME}/include)
set(ADIOS_LIBRARY ${ADIOS_HOME}/lib/libadios.a)
#set(MXML_LIBRARY /sw/redhat6/mxml/2.6/rhel6_gnu4.7.1/lib/libmxml.a)
set(MXML_LIBRARY /sw/redhat6/mxml/2.7/rhel6_gnu4.7.1/lib/libmxml.a)
set(ISOBAR_LIBRARY /ccs/proj/e2e/ncsu/sith.gnu/lib/libisobar.a)
set(APLOD_LIBRARY /ccs/proj/e2e/ncsu/sith.gnu/lib/libaplod.a) 
set(BZ2_LIBRARY /sw/sith/bzip2/1.0.6/rhel6_gnu4.7.1/lib/libbz2.a)
#set(DSPACES_LIBRARY ${DATASPACES_DIR}/lib/libdspaces.a)
#set(DSCOMMON_LIBRARY ${DATASPACES_DIR}/lib/libdscommon.a)
#set(DART_LIBRARY ${DATASPACES_DIR}/lib/libdart.a)
set(RDMACM_LIBRARY /usr/lib64/librdmacm.so)
#set(NETCDF_LIBRARY ${NETCDF_DIR}/lib/libnetcdf.a)
set(CURL_LIBRARY /usr/lib64/libcurl.so)
set(Z_LIBRARY /usr/lib64/libz.so)
set(M_LIBRARY /usr/lib64/libm.a) 
set(LUSTREAPI_LIBRARY /usr/lib64/liblustreapi.a)
set(IBVERBS_LIBRARY /usr/lib64/libibverbs.so)

set(CMAKE_C_COMPILER  gcc)
set(CMAKE_CXX_COMPILER  /sw/sith/ompi/1.6.5/rhel6_gnu4.7.1/bin/mpicxx)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER=1 -DUSE_REAL_STRUCT_FACTOR") 
set(GNU_FLAGS "-malign-double -fomit-frame-pointer -ffast-math -fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
#set(GNU_FLAGS "-g -O0 -malign-double -fomit-frame-pointer -ffast-math -fopenmp -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
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
  /sw/redhat6/hdf5-parallel/1.8.11/rhel6.4_gnu4.7.1
  /sw/redhat6/szip/2.1/rhel6_gnu4.7.1
  #/ccs/home/zgu/adioshub/adios_sith
  /ccs/proj/e2e/pnorbert/ADIOS/sith.gnu
  #/sw/redhat6/dataspaces/1.2.0/rhel6_gnu4.7.1
  /sw/sith/netcdf/4.1.3/rhel6_gnu4.7.1
  #/sw/redhat6/adios/1.6.0/rhel6_gnu4.7.1
  #/ccs/home/zgu/Software_sith
  /sw/redhat6/libxml2/2.9.1/rhel6_gnu4.7.1
  /sw/redhat6/fftw/3.3.3/rhel6.4_gnu4.7.1
  #/ccs/proj/e2e/pnorbert/qmcpack/Software_sith_cynthia/
)

link_libraries(/sw/redhat6/acml/5.3.0/gfortran64/lib/libacml.a -lz)
link_libraries(/ccs/compilers/gcc/rhel6-x86_64/4.7.1/lib64/libgfortran.so)
#link_libraries(-L/sw/redhat6/adios/1.6.0/rhel6_gnu4.7.1/lib -ladios -L/sw/sith/mxml/2.6/rhel6_gnu4.7.1/lib -lmxml -L/sw/redhat6/dataspaces/1.3.0/rhel6_gnu4.7.1/lib -L/usr/lib64/lib -L/lib64/lib -L/lib64/lib64 -L/lib -L/sw/redhat6/adios/1.6.0/rhel6_gnu4.7.1/lib -lm -lmxml -lpthread -ldspaces -ldscommon -ldart -lrdmacm -libverbs -llustreapi -lz -lbz2 -lisobar -lz)
#link_libraries(/ccs/home/zgu/adioshub/adios_sith/lib/libadios.a)
link_libraries(/sw/redhat6/libxml2/2.9.1/rhel6_gnu4.7.1/lib/libxml2.a)
link_libraries(/sw/redhat6/fftw/3.3.3/rhel6.4_gnu4.7.1/lib/libfftw3.a)
link_libraries(/ccs/proj/e2e/pnorbert/ADIOS/sith.gnu/lib/libadios.a)
link_libraries(/ccs/proj/e2e/ncsu/sith.gnu/lib/libisobar.a)
link_libraries(/ccs/proj/e2e/ncsu/sith.gnu/lib/libaplod.a)
link_libraries(/sw/sith/bzip2/1.0.6/rhel6_gnu4.7.1/lib/libbz2.a)
