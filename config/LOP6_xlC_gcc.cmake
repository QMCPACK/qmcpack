# tool chain for blueprint
SET(CMAKE_SYSTEM_PROCESSOR "P5p")
#SET(CMAKE_PLATFORM_REQUIRED_RUNTIME_PATH /usr/lib /lib)
SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG "-Wl,-blibpath:")
SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG_SEP ":")

# Files named "libfoo.a" may actually be shared libraries.
SET_PROPERTY(GLOBAL PROPERTY TARGET_ARCHIVES_MAY_BE_SHARED_LIBS 1)
#SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

set(CMAKE_C_COMPILER /sara/sw/modules11/wrappers/sara/gcc)
set(CMAKE_CXX_COMPILER /sara/sw/modules11/wrappers/sara/mpCC)
set(CMAKE_Fortran_COMPILER /sara/sw/modules11/wrappers/sara/xlf90)

set(CMAKE_FIND_ROOT_PATH
    /opt/ibmcmp/xlmass/5.0/
    /opt/ibmcmp/xlf/12.1/
    /home/jnkim/build/einspline/sles11_gnu
    /home/jnkim/build/fftw3.2/xlc64
    /home/jnkim/src/p3dfft/essl_stride1_20090928
    /sara/sw/hdf5/1.6.7
)

SET(AIX_ARCH_FLAGS "-g -q64 -qarch=auto -qtune=auto -qcache=auto -D_GNU_SOURCE -DUSE_EVEN -DIBM -DINLINE_ALL= ")
SET(AIX_OPT_FLAGS "-O3 -qmaxmem=-1 -qkeyword=restrict -qinline -qunroll=yes -qprefetch -qlargepage -qsmp=omp -qthreaded")
SET(GNU_ARCH_FLAGS "-g -mcpu=power6 -mtune=power6 -m64")
set(GNU_OPT_FLAGS "-fopenmp -O3 -Drestrict=__restrict__  -finline-functions -fexpensive-optimizations -ffast-math -finline-limit=1000 -fstrict-aliasing -funroll-all-loops")

#CXX flags
SET(AIX_CXX_FLAGS "-qsuppress=1540-1090:1540-1103:1540-1088:1540-0700")
#C flags
SET(AIX_C_FLAGS "-qlanglvl=stdc99")
set(GNU_C_FLAGS "-std=c99")

#INCLUDE(Platform/AIX)
######################################################################
#set the CXX flags: arch+common + opt 
######################################################################
SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS} ${AIX_CXX_FLAGS}")
SET(CMAKE_C_FLAGS   "${GNU_ARCH_FLAGS} ${GNU_OPT_FLAGS} ${GNU_C_FLAGS}")

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

SET(XLF_LIBS -L/opt/ibmcmp/xlf/12.1/lib64 -lxlf90_r)
link_libraries(-L/sara/sw/lapack/3.1.1/lib -llapack -lessl 
  -L/opt/ibmcmp/xlmass/5.0/lib64 -lmassvp6_64 -lmass_64  
  -L/opt/ibmcmp/vac/10.1/lib64 -lxl -lxlopt
  ${XLF_LIBS} -L/usr/lpp/mmfs/lib -lgpfs 
  -L/opt/ibmcmp/xlsmp/1.8/lib64 -lxlsmp 
  -lm)

SET(CMAKE_SKIP_RPATH TRUE)
