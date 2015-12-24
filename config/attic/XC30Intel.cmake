#######################################
# toolchain files on edison.nersc.gov
#######################################
# Issues 
# * Cannot get rid of -rdynamic for the linker and cmake fails badly.
#   need three steps (hack!!! hopefully some cmake guru can help me)
#  1. cmake -DCMAKE_TOOLCHAIN_FILE=~/svnwork/qmcpack-trunk/config/XC30.cmake ~/svnwork/qmcpack-trunk \
#     -DCMAKE_C_COMPILER_FORCED=true -DCMAKE_CXX_COMPILER_FORCED=true
#  2. cmake -DCMAKE_TOOLCHAIN_FILE=~/svnwork/qmcpack-trunk/config/XC30.cmake ~/svnwork/qmcpack-trunk
#  3. find . | grep link.txt | xargs -L1  sed -i s/\-rdynamic//
#  step 1 bypass compiler checks
#  step 3 remove -rdynamic
# * FMA is not working. Not critical, since float does not use it well anyway. 
# * BOOST is not picked up but workaround is below.
#
# So far
#   Other than the nightmare with cmake&-rdynamic, all is working well out of box. 
#
# Modules: add cmake & gcc (default gcc 4.3 and use 4.7.2)
#  module load cmake
#  module load gcc
#
# 4.3 may work but it is really old.
#
#Currently Loaded Modulefiles:
#  1) modules/3.2.6.6                      11) dmapp/4.0.1-1.0500.5932.6.5.ari      21) craype-sandybridge
#  2) nsg/1.2.0                            12) gni-headers/3.0-1.0500.5839.5.3.ari  22) cray-shmem/5.6.1
#  3) eswrap/1.0.19-1.010001.264.0         13) xpmem/0.1-2.0500.36799.3.6.ari       23) cray-mpich2/5.6.1
#  4) switch/1.0-1.0500.37046.2.39.ari     14) job/1.5.5-0.1_2.0500.35831.1.76.ari  24) torque/4.1.4
#  5) craype-network-aries                 15) csa/3.0.0-1_2.0500.37464.3.21.ari    25) moab/7.2.0-r1-b63-SUSE11
#  6) craype/1.01                          16) dvs/2.2_0.9.0-1.0500.1388.1.156.ari  26) altd/1.0
#  7) intel/13.0.1.117                     17) alps/5.0.1-2.0500.7663.1.1.ari       27) darshan/2.2.5-pre1
#  8) udreg/2.3.2-1.0500.5931.3.1.ari      18) rca/1.0.0-2.0500.37705.3.12.ari      28) usg-default-modules/1.0
#  9) ugni/4.0-1.0500.5836.7.58.ari        19) atp/1.6.1                            29) cmake/2.8.10.2
# 10) pmi/4.0.1-1.0000.9421.73.3.ari       20) PrgEnv-intel/5.0.15                  30) gcc/4.7.2

SET(CMAKE_SYSTEM_PROCESSOR "XC30")

set(CMAKE_CXX_COMPILER  /usr/common/usg/darshan/craype/1.06/bin/CC)
set(CMAKE_C_COMPILER  /usr/common/usg/darshan/craype/1.06/bin/cc)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline -DDISABLE_TIMER -DUSE_REAL_STRUCT_FACTOR")
set(INTEL_OPTS " -static -g  -restrict -unroll  -O3 -ip  -xAVX -openmp -Wno-deprecated")
set(XT_FLAGS "-D_CRAYMPI")# -DHAVE_FMA4=1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XT_FLAGS} ${INTEL_OPTS} ${GNU_OPTS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${XT_FLAGS} ${INTEL_OPTS} -std=c99")

set(CMAKE_C_COMPILER_FORCED)
#--------------------------------------------------------------------------
# path where the libraries are located
# boost,hdf,szip,libxml2,fftw,essl
#--------------------------------------------------------------------------
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_SHARED_LINKER_FLAGS) #-static -Wl,--export-dynamic")
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)# "-Wl,-Bstatic")
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)# "-Wl,-Bstatic")
set(CMAKE_EXE_LINKER_FLAGS "-static")

FOREACH(type SHARED_LIBRARY SHARED_MODULE EXE)
  #SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-static -Wl,--export-dynamic")
  #SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-static")
  #SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-static -Wl,--export-dynamic")
  #SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-static")
  SET(CMAKE_${type}_LINK_STATIC_C_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_C_FLAGS "-static")
  SET(CMAKE_${type}_LINK_STATIC_CXX_FLAGS "-Wl,-Bstatic")
  SET(CMAKE_${type}_LINK_DYNAMIC_CXX_FLAGS "-static")
ENDFOREACH(type)

set(CMAKE_FIND_ROOT_PATH
  /opt/cray/hdf5/1.8.9/intel/120
  /opt/fftw/3.3.0.1/x86_64
  /global/homes/j/jkim/share/xc30/libxml2
  /global/u2/j/jkim/share/boost_1_42_0
)

set(Boost_INCLUDE_DIR /global/u2/j/jkim/share/boost_1_42_0)
#--------------------------------------------------------------------------
# below is common for INTEL compilers and MKL library
#--------------------------------------------------------------------------
SET(QMC_BUILD_STATIC 1)
set(ENABLE_OPENMP 1)
set(HAVE_MPI 1)
set(HAVE_SSE 1)
set(HAVE_SSE2 1)
set(HAVE_SSE3 1)
set(HAVE_SSSE3 1)
set(HAVE_SSE41 1)
set(USE_PREFETCH 1)
set(PREFETCH_AHEAD 10)
set(HAVE_MKL 1)
set(HAVE_MKL_VML 1)

include_directories($ENV{MKLROOT}/include)
link_libraries(-L$ENV{MKLROOT}/lib/intel64 -mkl=sequential)

#link_libraries(-static -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm)
#link_libraries(/opt/intel/composer_xe_2013.1.117/mkl/lib/intel64/libmkl_sequential.a)
#INCLUDE(Platform/UnixPaths)
#
#SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
#SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
#SET(CMAKE_C_LINK_SHARED_LIBRARY)
#SET(CMAKE_C_LINK_MODULE_LIBRARY)
#
