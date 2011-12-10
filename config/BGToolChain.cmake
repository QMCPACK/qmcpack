# the name of the target operating system
#SET(CMAKE_SYSTEM_NAME BlueGeneP)
SET(BGP 1 CACHE BOOL "On BlueGeneP")
SET(Linux 0)

SET(QMC_BUILD_STATIC 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_LIBESSL 1)

# set the compiler
set(CMAKE_C_COMPILER  /soft/apps/darshan/bin/default/mpixlc_r)
set(CMAKE_CXX_COMPILER  /soft/apps/darshan/bin/default/mpixlcxx_r)

set(AIX_ARCH "450d")
SET(AIX_ARCH_FLAGS "-qarch=${AIX_ARCH}  -qsmp=omp:noopt -qthreaded -qkeyword=restrict -qstrict")
#SET(AIX_CXX_COMMON_FLAGS " -qkeyword=restrict -qstrict -qhot -qnoeh -qsuppress=1540-1090:1540-1088 ")
#SET(AIX_CXX_COMMON_FLAGS " -qkeyword=restrict -qstrict -qnoeh -qsuppress=1540-1090:1540-1088 ")
SET(AIX_CXX_COMMON_FLAGS "-qdebug=NADDRTKNCFG -qsuppress=1540-1090:1540-1088 -qpath=ILbc:/soft/apps/xlC-interim-fix-ICE-qmcpack/exe/")
SET(AIX_OPT_FLAGS "-qmaxmem=-1 -qprefetch ")

#SET(AIX_OPT_FLAGS "-O3 -Q -qmaxmem=-1 -qipa=inline -qinline -qlargepage -qprefetch ")
#SET(AIX_CXX_OPT_FLAGS "-O3 -Q -qlargepage -qprefetch")
#SET(AIX_CXX_FLAGS "-O3 -Q -qlargepage -qprefetch")

SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS}")

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")
ENDFOREACH(type)

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /home/projects/qmcpack/boost-1.45_0
    /home/projects/qmcpack/einspline
    /home/projects/qmcpack/libxml2
    /soft/apps/hdf5-1.6.6
    /soft/apps/fftw-3.1.2-double
)

set(LAPACK_LIBRARY /soft/apps/LAPACK/lapack_3.3_BGP.a)
set(BLAS_LIBRARY /soft/apps/ESSL-4.4.1-1/lib/libesslbg.a)
set(FORTRAN_LIBRARIES 
/soft/apps/ibmcmp-aug2011/xlf/bg/11.1/bglib/libxlf90_r.a
/soft/apps/ibmcmp-aug2011/xlf/bg/11.1/bglib/libxlfmath.a
/soft/apps/ibmcmp-aug2011/xlf/bg/11.1/bglib/libxlopt.a
)
#link lapack, essl, mass
#  link_libraries(/soft/apps/LAPACK/lapack_3.3_BGP.a 
#      /soft/apps/ESSL/lib/libesslbg.a
#      /soft/apps/ibmcmp-aug2011/xlf/bg/11.1/lib/libxlf90_r.a 
#      /soft/apps/ibmcmp-aug2011/xlf/bg/11.1/lib/libxlfmath.a
#      /soft/apps/ibmcmp-aug2011/xlf/bg/11.1/lib/libxlopt.a
#      )
#
