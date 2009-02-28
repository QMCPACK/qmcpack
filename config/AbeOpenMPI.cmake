# tool chain for abe Fri Feb 27 2009
SET(CMAKE_SYSTEM_PROCESSOR "ES")
SET(QMC_ENV "LinuxIntel" CACHE STRING "Setting envirnoments for abe@ncsa")
SET(MKL_HOME "/usr/local/intel/mkl/10.0.3.020" CACHE STRING "MKL HOME")
SET(MKL_EXT "em64t" CACHE STRING "MKL Extension")
SET(EINSPLINE_HOME "/u/ac/esler" CACHE STRING "Einspline library directory")
#SET(EINSPLINE_HOME "/u/ac/jnkim/svnwork/einspline" CACHE STRING "Einspline source directory")

set(CMAKE_C_COMPILER  /usr/local/intel/10.1.017/bin/icc)
set(CMAKE_CXX_COMPILER /usr/local/openmpi-1.2.4-intel-ofed-1.2/bin/mpicxx)

set(CMAKE_FIND_ROOT_PATH
    /u/ac/jnkim/share/boost
    /usr/apps/hdf/hdf5/v167
    /usr/apps/hdf/szip/v2.1/static/encoder
    /usr/local/libxml2-2.6.29
    /usr/apps/math/fftw/fftw-3.1.2/intel10
    /u/ac/esler
    )

set(INTEL_OPTS "-g  -restrict -unroll  -O3 -ip -xT -openmp -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_OPTS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_OPTS} -std=c99")

INCLUDE(Platform/UnixPaths)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)
