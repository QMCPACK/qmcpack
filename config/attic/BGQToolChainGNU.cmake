# the name of the target operating system
#SET(CMAKE_SYSTEM_NAME BlueGeneP)
SET(BGP 0 CACHE BOOL "On BlueGeneP")
SET(Linux 0)

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_LIBESSL 1)
SET(HAVE_EINSPLINE 1)
SET(HAVE_EINSPLINE_EXT 0)
set(HAVE_ADIOS 0)




# set the compiler
set(CMAKE_C_COMPILER mpicc) 
set(CMAKE_CXX_COMPILER mpicxx) 

set(GNU_OPTS " -DINLINE_ALL=inline")
set(GNU_FLAGS "-fopenmp -g -O3 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -static -Wno-deprecated -Wl,-allow-multiple-definition ")

SET(CMAKE_CXX_FLAGS " ${GNU_FLAGS} ${GNU_OPTS} -DHAVE_MASSV -DHAVE_MASS -DSPLINEFLOAT")
#SET(CMAKE_CXX_FLAGS " ${GNU_FLAGS} ${GNU_OPTS} -DHAVE_MASSV -DHAVE_MASS ")
SET(CMAKE_C_FLAGS " -DHAVE_MASSV -DHAVE_MASS -std=gnu99  -fomit-frame-pointer ${GNU_FLAGS} ${GNU_OPTS} ")

## use this to search only these directores not standard directories.
#include_directories(/bgsys/drivers/ppcfloor/comm/gcc/include)
set(Boost_INCLUDE_DIR /home/projects/qmcpack/boost_1_45_0)
SET(CMAKE_FIND_ROOT_PATH  
     /home/projects/qmcpack/boost_1_45_0
     /home/projects/qmcpack/LIBXML2-2.9
     /soft/libraries/hdf5/1.8.10/cnk-xl/current
     /soft/libraries/alcf/current/xl/ZLIB
     /soft/libraries/alcf/current/gcc/FFTW3
     /soft/compilers/ibmcmp-nov2013/xlmass/bg/7.3/include 


)
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)


include_directories(/soft/compilers/ibmcmp-nov2013/xlmass/bg/7.3/include 

)
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)


include_directories(/soft/compilers/ibmcmp-nov2013/xlmass/bg/7.3/include
#/soft/perftools/bgpm/include/
)


set(LAPACK_LIBRARY /soft/libraries/alcf/current/xl/LAPACK/lib/liblapack.a)  
set(BLAS_LIBRARY /soft/libraries/essl/5.1.1-0.beta/lib64/libesslbg.a)
SET(FORTRAN_LIBRARIES 
/soft/compilers/ibmcmp-nov2013/xlf/bg/14.1/bglib64/libxlf90_r.a
/soft/compilers/ibmcmp-nov2013/xlf/bg/14.1/bglib64/libxlfmath.a
/soft/compilers/ibmcmp-nov2013/xlf/bg/14.1/bglib64/libxlopt.a
/soft/compilers/ibmcmp-nov2013/xlf/bg/14.1/bglib64/libxl.a
/bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux/lib/libdl.a
)
link_libraries(/bgsys/drivers/V1R2M0/ppc64/gnu-linux/powerpc64-bgq-linux/lib/libgfortran.a
/soft/compilers/ibmcmp-nov2013/xlmass/bg/7.3/bglib64/libmass.a 
/soft/compilers/ibmcmp-nov2013/xlmass/bg/7.3/bglib64/libmassv.a
#/soft/perftools/hpctw/libmpihpm_smp.a
#/bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a
#/bgsys/drivers/ppcfloor/spi/lib/libSPI_upci_cnk.a
#-pg
)

SET(BUILD_SANDBOX 1)


