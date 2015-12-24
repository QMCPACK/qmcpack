# tool chain for blueprint
SET(CMAKE_SYSTEM_PROCESSOR "P5p")
#SET(CMAKE_PLATFORM_REQUIRED_RUNTIME_PATH /usr/lib /lib)
SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG "-Wl,-blibpath:")
SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG_SEP ":")

# Files named "libfoo.a" may actually be shared libraries.
SET_PROPERTY(GLOBAL PROPERTY TARGET_ARCHIVES_MAY_BE_SHARED_LIBS 1)
#SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

set(CMAKE_C_COMPILER /sara/sw/modules11/wrappers/sara/mpcc)
set(CMAKE_CXX_COMPILER /sara/sw/modules11/wrappers/sara/mpCC)
set(CMAKE_Fortran_COMPILER /sara/sw/modules11/wrappers/sara/xlf90)

set(CMAKE_FIND_ROOT_PATH
    /opt/ibmcmp/xlmass/5.0/
    /opt/ibmcmp/xlf/12.1/
    /home/jnkim/build/fftw3.2/xlc64
    /home/jnkim/src/p3dfft/essl_stride1_20090928
    /opt/ibmhpc/ppe.hpct
    /sara/sw/hdf5/1.8.2pp
     )

#SET(AIX_ARCH_FLAGS "-g -q64 -qarch=auto -qtune=auto -qcache=auto -DUSE_EVEN -DIBM -DUSE_ALLTOALLV")
SET(AIX_ARCH_FLAGS "-g -q64 -qarch=auto -qtune=auto -qcache=auto -DUSE_EVEN -DIBM -DINLINE_ALL= ")
SET(AIX_CXX_COMMON_FLAGS "-qnoeh -qsuppress=1540-1090:1540-1103:1540-1088:1540-0700 -qsmp=omp")
#SET(AIX_OPT_FLAGS "-O3 -Q -qmaxmem=-1 -qinline -qlargepage -qprefetch -qstrict -qhot -qkeyword=restrict")
#SET(AIX_F_FLAGS "-O3 -Q -qmaxmem=-1 -qinline -qlargepage -qprefetch -qstrict -qhot")
SET(AIX_OPT_FLAGS "-O3 -qmaxmem=-1 -qprefetch -qstrict -qhot -qkeyword=restrict")
SET(AIX_F_FLAGS "-O3 -qmaxmem=-1 -qprefetch -qstrict -qhot")

#SET(AIX_OPT_FLAGS "-O3 -Q -qmaxmem=-1 -qipa=inline -qinline -qlargepage -qprefetch -qstrict -qhot -qkeyword=restrict")
#SET(AIX_OPT_FLAGS "-g -O3 -qlargepage -qprefetch -qstrict -qhot -qkeyword=restrict")
#SET(AIX_OPT_FLAGS "-g -Q -qmaxmem=-1 -qlargepage -qprefetch -qkeyword=restrict")
#SET(AIX_C_FLAGS "-qlanglvl=stdc99")

#INCLUDE(Platform/AIX)
######################################################################
#set the CXX flags: arch+common + opt 
######################################################################
SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS} ${AIX_C_FLAGS}")
SET(CMAKE_Fortran_FLAGS "${F_DEFINES} ${AIX_ARCH_FLAGS} ${AIX_F_FLAGS}")
SET(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})

SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_LIBESSL 1)
SET(ENABLE_FORTRAN 1)

SET(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "-G -Wl,-brtl,-bnoipath")  # -shared
SET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-Wl,-brtl,-bnoipath,-bexpall")  # +s, flag for exe link to use shared lib
SET(CMAKE_SHARED_LIBRARY_Fortran_FLAGS " ")
SET(CMAKE_SHARED_MODULE_Fortran_FLAGS  " ")

# CXX Compiler
IF(CMAKE_COMPILER_IS_GNUCXX) 
  SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "-shared -Wl,-G")       # -shared
ENDIF(CMAKE_COMPILER_IS_GNUCXX) 

# C Compiler
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "-shared -Wl,-G")       # -shared
ENDIF(CMAKE_COMPILER_IS_GNUCC)

SET(XLF_LIBS  -L/opt/ibmcmp/xlf/12.1/lib64 -lxlf90_r -lxlomp_ser -lxlsmp -lxl -lxlopt)
link_libraries(-llapack -lessl -lmassvp6_64 -lmass_64  ${XLF_LIBS} -L/usr/lpp/mmfs/lib -lgpfs)
#SET(BLAS_LIBRARY -lessl -lmassvp6_64 -lmass_64 ${XLF_LIBS})
#SET(LAPACK_LIBRARY -llapack)

SET(CMAKE_SKIP_RPATH TRUE)
