SET(CMAKE_SYSTEM_PROCESSOR "x86_64")
# Last updated: July 2, 2015
# The following modules are required
# Currently Loaded Modulefiles:
#   1) libxml/2.9.2                         5) /numactl/2.0.9
#   2) gcc/4.8.4                            6) mvapich2-2.1rc2/gcc-4.8.4/cuda-7.0
#   3) /fftw/3.3.3-withthreads              7) intel/mkl/11.0
#   4) cuda/7.0
# The libxml module can be found after the following command
# module use /home-2/jlarkin/sw/modulefiles
#
# Build inside an interactive PBS job to avoid "invalid instruction" errors. 
# This has been tested using the following PBS submission command.
#
# qsub -lnodes=1:ppn=20:k40:hsw

set(CMAKE_C_COMPILER  mpicc)
set(CMAKE_CXX_COMPILER  mpic++)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS " -fomit-frame-pointer -malign-double -fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops ")
set(OPT_FLAGS " -msse -msse2 -msse3 -msse4.1")
set(CMAKE_CXX_FLAGS "${OPT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS} -Wno-deprecated")
set(CMAKE_C_FLAGS "${OPT_FLAGS} ${GNU_FLAGS} -std=c99")

# need for both c++ and c
set(QMC_CUDA 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
set(HAVE_CUDA 1)
SET(QMC_BUILD_STATIC 1)
message(STATUS "QMC_CUDA: " ${QMC_CUDA})

link_libraries(-L/shared/apps/rhel-6.2/intel/ics-2013/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential
  -lmkl_core)

# Are these the correct settings for a Haswell CPU?
SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
SET(HAVE_SSE41 1)
SET(USE_PREFETCH 1)
SET(PREFETCH_AHEAD 12)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_SHARED_LINKER_FLAGS "")

set(CMAKE_FIND_ROOT_PATH
  /shared/apps/rhel-6.2/libs/hdf5-1.8.9/serial/
  /shared/apps/rhel-6.2/intel/ics-2013/mkl/lib/intel64
  $ENV{FFTW_DIR}/..
  $ENV{BOOST_DIR}
)

# bypass einspline search
#set(EINSPLINE_HOME /lustre/widow3/scratch/jnkim/einspline)
#set(HAVE_EINSPLINE 1)
#set(HAVE_EINSPLINE_EXT 0)

set(CUDA_NVCC_FLAGS "-Xcompiler=-fpermissive;-arch=sm_35;-Drestrict=__restrict__;-DNO_CUDA_MAIN;-O3")
set(CUDA_CUDART_LIBRARY "-lcuda")
set(CUDA_CUDA_LIBRARY "-cuda")
set(CUDA_cublas_LIBRARY $ENV{CUDA_HOME}/lib64/libcublas.so)
set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_HOME})
set(CUDA_TOOLKIT_INCLUDE $ENV{CUDA_HOME}/include)
set(CUDA_LIBRARIES ${CUDA_CUDART_LIBRARY})
set(CUDA_make2cmake ${CMAKE_ROOT}/Modules/FindCUDA/make2cmake.cmake)
set(CUDA_parse_cubin ${CMAKE_ROOT}/Modules/FindCUDA/parse_cubin.cmake)
set(CUDA_run_nvcc ${CMAKE_ROOT}/Modules/FindCUDA/run_nvcc.cmake)
set(CUDA_PROPAGATE_HOST_FLAGS OFF) #do not propagate the host flags


message(STATUS "CUDA_LIBRARY" $ENV{CUDA_CUDART_LIBRARY})

