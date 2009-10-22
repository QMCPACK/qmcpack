# tool chain for blueprint
SET(CMAKE_SYSTEM_PROCESSOR "P6")
#SET(QMC_ENV "IBMP5p" CACHE STRING "Setting envirnoments for IBM P5+")

SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

set(CMAKE_C_COMPILER /opt/ibmcmp/vacpp/10.1/bin/xlc_r)
set(CMAKE_CXX_COMPILER /opt/ibmhpc/ppe.poe/bin/mpCC)
#set(CMAKE_CXX_COMPILER /sara/sw/modules/wrappers/sara/xlC_r)
set(CMAKE_Fortran_COMPILER /sara/sw/modules/wrappers/sara/xlf_r)

set(CMAKE_FIND_ROOT_PATH
    /opt/ibmcmp/xlmass/5.0/
    /opt/ibmcmp/xlf/12.1/
    /sara/sw/hdf5/1.8.3
    /home/jnkim/build/fftw3.2/xlc64
     )

#SET(AIX_ARCH_FLAGS "-g -q64 -qarch=auto -qtune=auto -qcache=auto -DUSE_EVEN -DIBM -DUSE_ALLTOALLV")
SET(AIX_ARCH_FLAGS "-g -q64 -qarch=auto -qtune=auto -qcache=auto -DIBM -DHAVE_MASSV -DINLINE_ALL= ")
SET(AIX_CXX_COMMON_FLAGS "-qnoeh -qsuppress=1540-1090:1540-1103:1540-1088 -qsmp=omp")
#SET(AIX_OPT_FLAGS "-O3 -Q -qmaxmem=-1 -qinline -qlargepage -qprefetch -qstrict -qhot -qkeyword=restrict")
#SET(AIX_F_FLAGS "-O3 -Q -qmaxmem=-1 -qinline -qlargepage -qprefetch -qstrict -qhot")
SET(AIX_OPT_FLAGS "-O3 -qmaxmem=-1 -qprefetch -qlargepage -qhot -qkeyword=restrict")
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

SET(XLF_LIBS  -L/opt/ibmcmp/xlf/12.1/lib64 -lxlf90_r -lxl -lxlopt)
#SET(BLAS_LIBRARY -lessl -lmassvp6_64 -lmass_64 -L/usr/local/ihpct_2.2/lib64 -lhpm ${XLF_LIBS})
#SET(BLAS_LIBRARY -lessl -lmassvp6_64 -lmass_64 -L/usr/local/ihpct_2.2/lib64 -lmpitrace -llicense ${XLF_LIBS})
SET(BLAS_LIBRARY -lessl -lmassvp6_64 -lmass_64 ${XLF_LIBS})
SET(LAPACK_LIBRARY -L/sara/sw/lapack/3.1.1/lib -llapack)


#SET(FORTRAN_LIBS  -lxlf90 -lxlf -lxlopt)
#SET(CMAKE_CXX_CREATE_STATIC_LIBRARY
#    "<CMAKE_AR> -X64 cr <TARGET> <LINK_FLAGS> <OBJECTS> " 
#    "<CMAKE_RANLIB> <TARGET> ")
###SET(CMAKE_SHARED_LIBRARY_PREFIX "lib")          # lib
##SET(CMAKE_SHARED_LIBRARY_SUFFIX ".so")          # .so
##SET(CMAKE_DL_LIBS "-lld")
##
### RPATH support on AIX is called libpath.  By default the runtime
### libpath is paths specified by -L followed by /usr/lib and /lib.  In
### order to prevent the -L paths from being used we must force use of
### -Wl,-blibpath:/usr/lib:/lib whether RPATH support is on or not.
### When our own RPATH is to be added it may be inserted before the
### "always" paths.
##SET(CMAKE_PLATFORM_REQUIRED_RUNTIME_PATH /usr/lib /lib)
##SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG "-Wl,-blibpath:")
##SET(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG_SEP ":")
##
### Files named "libfoo.a" may actually be shared libraries.
##SET_PROPERTY(GLOBAL PROPERTY TARGET_ARCHIVES_MAY_BE_SHARED_LIBS 1)
##
### CXX Compiler
IF(CMAKE_COMPILER_IS_GNUCXX) 
  SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "-shared -Wl,-G")       # -shared
ELSE(CMAKE_COMPILER_IS_GNUCXX) 
  SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS " ")       # -shared
ENDIF(CMAKE_COMPILER_IS_GNUCXX) 
#  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Wl,-brtl,-bnoipath,-bexpall")  # +s, flag for exe link to use shared lib
#Safe choice
  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-B/usr/share/libhugetlbfs/ -tl -Wl,--hugetlbfs-link=BDT -L/usr/lib64 ")  # +s, flag for exe link to use shared lib
#link hugh page
#  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-B/usr/share/libhugetlbfs/ -tl -Wl,--hugetlbfs-link=BDT  -L/usr/lib64")  # +s, flag for exe link to use shared lib
  SET(CMAKE_SHARED_LIBRARY_CXX_FLAGS " ")
  SET(CMAKE_SHARED_MODULE_CXX_FLAGS  " ")
  SET (CMAKE_CXX_FLAGS_DEBUG_INIT "-g")
  SET (CMAKE_CXX_FLAGS_RELEASE_INIT "-O -DNDEBUG") 
  SET (CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O -DNDEBUG")
  SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-g")
