# tool chain for blueprint
SET(CMAKE_SYSTEM_PROCESSOR "P7")
SET(QMC_ENV "GCC4.5" CACHE STRING "Setting envirnoments for IBM P7")

SET(HAVE_MPI 0)
SET(ENABLE_OPENMP 1)
SET(HAVE_ESSL 1)
SET(IBM_COMPILER 0)
SET(HAVE_LIBHDF5 1)
SET(HAVE_EINSPLINE 1)
SET(HAVE_EINSPLINE_EXT 1)

set(CMAKE_FIND_ROOT_PATH
    /home/jnkim/share/boost
    /opt/hdf5-1.8.4-patch1-xl64
    /home/jnkim/share/at40/einspline
    /opt/fftw-3.2.2
)

#SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)
add_definitions(-DINLINEALL=inline -DINLINE0=inline -DINLINE1=inline )
set(CMAKE_C_COMPILER /opt/ibmcmp/vacpp/11.1/bin/xlc_r)
set(CMAKE_CXX_COMPILER /opt/at4.0/bin/g++)
set(CMAKE_Fortran_COMPILER /opt/ibmcmp/xlf/13.1/bin/xlf_r)

SET(AIX_ARCH_FLAGS "-q64 -qarch=pwr7 -qtune=pwr7 -qcache=auto -qsimd=auto -qhot=simd -qinline -qsmp=omp -qthreaded -unroll=yes")
SET(AIX_OPT_FLAGS "-O3 -qmaxmem=-1 -qlargepage -qprefetch -qkeyword=restrict")
SET(F_DEFINES "-WF,-P,-DESSL,-DDOUBLE_PREC,-DSTRIDE1 -qfixed=132")
SET(AIX_F_FLAGS "-O3 -Q -qmaxmem=-1 -qlargepage -qprefetch")

SET(GNU_ARCH_FLAGS "-g -m64 -mcpu=power7 -mtune=power7 -mvsx -DADD_")
set(GNU_OPT_FLAGS "-fopenmp -O3 -Drestrict=__restrict__  -finline-functions -fexpensive-optimizations -ffast-math -finline-limit=1000 -fstrict-aliasing -funroll-all-loops")
set(GNU_CXX_FLAGS " -ftemplate-depth-60 ")

######################################################################
#set the CXX flags: arch+common + opt 
######################################################################
#SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
SET(CMAKE_CXX_FLAGS "${GNU_ARCH_FLAGS} ${GNU_OPT_FLAGS} ${GNU_CXX_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS} ${AIX_C_FLAGS}")
SET(CMAKE_Fortran_FLAGS "${F_DEFINES} ${AIX_ARCH_FLAGS} ${AIX_F_FLAGS}")
SET(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})

SET(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "-G -Wl,-brtl,-bnoipath")  # -shared
SET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-Wl,-brtl,-bnoipath,-bexpall")  # +s, flag for exe link to use shared lib
SET(CMAKE_SHARED_LIBRARY_Fortran_FLAGS " ")
SET(CMAKE_SHARED_MODULE_Fortran_FLAGS  " ")
#SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Wl,-brtl,-bnoipath,-bexpall")  # +s, flag for exe link to use shared lib
#SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Wl,--start-group")

# CXX Compiler
IF(CMAKE_COMPILER_IS_GNUCXX) 
  SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "-shared -Wl,-G")       # -shared
ENDIF(CMAKE_COMPILER_IS_GNUCXX) 

# C Compiler
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "-shared -Wl,-G")       # -shared
ENDIF(CMAKE_COMPILER_IS_GNUCC)

INCLUDE_DIRECTORIES(/opt/ibmcmp/xlmass/6.1/include)

link_libraries(
  -Wl,--start-group
  -L/home/jnkim/share/at40/lib -llapack  -lblas -lgfortran
  -L/opt/ibmcmp/xlmass/6.1/lib64 -lmass_simdp7_64 -lmassvp7_64
  /opt/hdf5-1.8.4-patch1-xl64/lib/libhdf5.a
  -Wl,--end-group
  )
