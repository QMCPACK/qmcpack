# tool chain for blueprint
SET(CMAKE_SYSTEM_PROCESSOR "P7")
SET(QMC_ENV "IBMP7" CACHE STRING "Setting envirnoments for IBM P5+")

SET(HAVE_MPI 1)
SET(ENABLE_OPENMP 1)
SET(ENABLE_FORTRAN 0)
SET(HAVE_ESSL 1)


set(CMAKE_C_COMPILER  /opt/ibmcmp/vacpp/11.1/bin/xlc_r)
#set(CMAKE_CXX_COMPILER /opt/ibmcmp/vacpp/11.1/bin/xlC_r)
set(CMAKE_CXX_COMPILER /usr/bin/mpCC)
set(CMAKE_Fortran_COMPILER /opt/ibmcmp/xlf/13.1/bin/xlf_r)

SET(AIX_ARCH_FLAGS "-q64 -qarch=pwr7 -qtune=pwr7 -qcache=auto -qsimd=auto -qhot=simd -qsmp=omp -qthreaded")
SET(AIX_CXX_COMMON_FLAGS "-qsuppress=1540-1090:1540-1088 -DINLINE0= -DINLINE1= -DDOUBLE_PREC")
SET(AIX_OPT_FLAGS "-O3 -qstrict -qkeyword=restrict -qmaxmem=-1 -qunroll=yes -qprefetch")
#SET(AIX_OPT_FLAGS "-O3 -qhot -qkeyword=restrict -qmaxmem=-1 -qlargepage -qunroll=yes -qprefetch")

#SET(AIX_OPT_FLAGS "-O3 -Q -qmaxmem=-1 -qipa=inline -qinline -qlargepage -qprefetch -Drestrict=__restrict__ -qkeyword=restrict -qunroll=yes")
#SET(AIX_OPT_FLAGS "-O3 -qmaxmem=-1 -qlargepage -qprefetch -Drestrict=__restrict__ -qkeyword=restrict -qunroll=yes")
#SET(AIX_OPT_FLAGS "-O3 -qhot=level=1 -qkeyword=restrict -qmaxmem=-1 -qlargepage -qprefetch=assistthread -qunroll=yes ")
#SET(AIX_OPT_FLAGS "-O3 -qhot -qkeyword=restrict -qmaxmem=-1 -qlargepage -qprefetch=assistthread -qunroll=yes ")
#SET(AIX_OPT_FLAGS "-O3 -qinline -qhot -qkeyword=restrict -qmaxmem=-1 -qlargepage -qunroll=yes -qprefetch")
#SET(F_DEFINES "-WF,-P,-DESSL,-DDOUBLE_PREC,-DSTRIDE1 -qfixed=132")
#SET(AIX_OPT_FLAGS "-g -O3 -qlargepage -qprefetch -qstrict -qhot -qkeyword=restrict")
#SET(AIX_OPT_FLAGS "-g -Q -qmaxmem=-1 -qlargepage -qprefetch -qkeyword=restrict")
#SET(AIX_C_FLAGS "-qlanglvl=stdc99")


SET(F_DEFINES "-qfixed=132")
SET(AIX_F_FLAGS "-O3 -qmaxmem=-1 -qlargepage -qhot -qunroll=yes")

set(CMAKE_FIND_ROOT_PATH
  /u02/ncsa03/share/xl/libxml
  /u02/ncsa03/share/xl/hdf5
  /u02/ncsa03/share/xl/fftw-3.2.2
  /u02/ncsa03/share/xl/einspline
  /u02/ncsa03/share/boost
  )

######################################################################
#set the CXX flags: arch+common + opt 
######################################################################
SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS} ${AIX_C_FLAGS}")
SET(CMAKE_Fortran_FLAGS "${F_DEFINES} ${AIX_ARCH_FLAGS} ${AIX_F_FLAGS}")
SET(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS})
SET(MAKE_Fortran_IMPLICIT_LINK_LIBRARIES "xlf90_r")
#SET(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortranbegin;gfortran;...")
#SET(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "...") 

#using LAPACK/ESSL
link_libraries(-L/u02/ncsa04/lib64 -llapack -L/usr/lib64 -lessl -L/opt/ibmcmp/xlmass/6.1/lib64  -lmass_simdp7_64 -L/opt/ibmcmp/xlf/13.1/lib64 -lxlf90_r -lxlsmp)
#using LAPACK/BLAS
#link_libraries(-L/u02/ncsa04/lib64 -llapack -lblas -L/opt/ibmcmp/xlmass/6.1/lib64  -lmass_simdp7_64 -L/opt/ibmcmp/xlf/13.1/lib64 -lxlf90_r -lxlsmp)


#link_libraries(-lessl -lmass -lmassv -lxlf90_r)
#link_libraries(-lessl -lmass -lmassv -lxlf90_r  -L/usr/lpp/ppe.hpct/lib64 -lmpitrace)
#link_libraries(-lessl -lmass -lmassv -lxlf90_r -L/u02/ncsa04/lib -lmpiP)
#link_libraries(-L/usr/lib64 -lessl  -L/u02/ncsa04/lib64 -llapack -L/opt/ibmcmp/xlmass/6.1/lib64  -lmass_simdp7_64 -L/opt/ibmcmp/xlf/13.1/lib64 -lxlf90_r -lxlsmp)
#link_libraries(-L/opt/apps/ibm/essl/5.1/lib64 -lessl -L/opt/apps/ibm/xlf/13.1/lib64 -lxlf90_r)
