# tool chain for blueprint
SET(CMAKE_SYSTEM_PROCESSOR "P5p")

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
    /sara/sw/hdf5/1.8.3
     )

   #set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_OPTS " -DINLINE_ALL=inline")
set(GNU_FLAGS "-fopenmp -O3 -Drestrict=__restrict__  -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
#SET(AIX_ARCH_FLAGS "-g -q64 -qarch=auto -qtune=auto -qcache=auto -DUSE_EVEN -DIBM -DUSE_ALLTOALLV")
SET(AIX_ARCH_FLAGS "-m64 -g")
set(AIX_CXX_FLAGS "-compiler g++ -ftemplate-depth-60")

SET(XLC_FLAGS "-g -q64 -qarch=auto -qtune=auto -qcache=auto -O3 -qmaxmem=-1 -qprefetch -qstrict -qhot -qkeyword=restrict")

#INCLUDE(Platform/AIX)
######################################################################
#set the CXX flags: arch+common + opt 
######################################################################
SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${GNU_FLAGS} ${GNU_OPTS} ${AIX_CXX_FLAGS}")
SET(CMAKE_C_FLAGS   "${XLC_FLAGS}")
#SET(CMAKE_Fortran_FLAGS "${F_DEFINES} ${AIX_ARCH_FLAGS} ${AIX_F_FLAGS}")
#SET(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})

SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
SET(HAVE_LIBESSL 1)
SET(ENABLE_FORTRAN 1)

#SET(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "-G -Wl,-brtl,-bnoipath")  # -shared
#SET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-Wl,-brtl,-bnoipath,-bexpall")  # +s, flag for exe link to use shared lib
#SET(CMAKE_SHARED_LIBRARY_Fortran_FLAGS " ")
#SET(CMAKE_SHARED_MODULE_Fortran_FLAGS  " ")

# CXX Compiler
IF(CMAKE_COMPILER_IS_GNUCXX) 
  SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "-shared -Wl,-G")       # -shared
ENDIF(CMAKE_COMPILER_IS_GNUCXX) 

# C Compiler
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "-shared -Wl,-G")       # -shared
ENDIF(CMAKE_COMPILER_IS_GNUCC)

#SET(XLF_LIBS  -L/opt/ibmcmp/xlf/12.1/lib64 -lxlf90_r -lxlomp_ser -lxlsmp -lxl -lxlopt)
SET(XLC_LIBS -L/sara/sw/lapack/3.1.1/lib -llapack -lessl
  -L/opt/ibmcmp/vac/10.1/lib64 -lxl -lxlopt 
  -L/opt/ibmcmp/xlsmp/1.8/lib64 -lxlsmp -lxlomp_ser 
  -L/opt/ibmcmp/xlf/12.1/lib64 -lxlf90_r 
  -L/opt/ibmcmp/xlmass/5.0/lib64 -lmassvp6_64 -lmass_64 
  -L/usr/lpp/mmfs/lib -lgpfs)

#SET(CMAKE_SKIP_RPATH TRUE)
