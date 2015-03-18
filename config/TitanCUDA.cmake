SET(CMAKE_SYSTEM_PROCESSOR "XK7")
# Last updated: Nov 7, 2014
# NEED THESES + defaults
#  module swap PrgEnv-pgi PrgEnv-gnu
#  module load cray-hdf5
#  module load fftw
#  module load boost
#  module load subversion
#  module load cmake
#  module load cudatoolkit
set(CMAKE_C_COMPILER  cc)
set(CMAKE_CXX_COMPILER  CC)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS " -fomit-frame-pointer -malign-double -fopenmp -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops ")
#set(XT_FLAGS "-march=amdfam10 -msse3 -D_CRAYMPI")
#set(XT_FLAGS "-march=bdver1 -msse3 -D_CRAYMPI")
set(XT_FLAGS " -msse -msse2 -msse3 -msse4.1 -D_CRAYMPI")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS} -Wno-deprecated")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -std=c99")

# need for both c++ and c
set(QMC_CUDA 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
set(HAVE_CUDA 1)
SET(QMC_BUILD_STATIC 1)
message(STATUS "QMC_CUDA: " ${QMC_CUDA})

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
  $ENV{HDF5_DIR}
  $ENV{FFTW_DIR}/..
  $ENV{BOOST_DIR}
  /sw/xk7/libxml2
)

# bypass einspline search
#set(EINSPLINE_HOME /lustre/widow3/scratch/jnkim/einspline)
#set(HAVE_EINSPLINE 1)
#set(HAVE_EINSPLINE_EXT 0)

set(CUDA_NVCC_FLAGS "-arch=sm_35;-Drestrict=__restrict__;-DNO_CUDA_MAIN;-O3")
set(CUDA_CUDART_LIBRARY /opt/cray/nvidia/default/lib64/libcuda.so)
set(CUDA_CUDA_LIBRARY /opt/cray/nvidia/default/lib64/libcuda.so)
set(CUDA_cublas_LIBRARY $ENV{CRAY_CUDATOOLKIT_DIR}/lib64/libcublas.so)
set(CUDA_TOOLKIT_ROOT_DIR $ENV{CRAY_CUDATOOLKIT_DIR})
set(CUDA_TOOLKIT_INCLUDE $ENV{CRAY_CUDATOOLKIT_DIR}/include)
set(CUDA_LIBRARIES ${CUDA_CUDART_LIBRARY})
set(CUDA_make2cmake ${CMAKE_ROOT}/Modules/FindCUDA/make2cmake.cmake)
set(CUDA_parse_cubin ${CMAKE_ROOT}/Modules/FindCUDA/parse_cubin.cmake)
set(CUDA_run_nvcc ${CMAKE_ROOT}/Modules/FindCUDA/run_nvcc.cmake)
set(CUDA_PROPAGATE_HOST_FLAGS OFF) #do not propagate the host flags


message(STATUS "CUDA_LIBRARY" $ENV{CUDA_CUDART_LIBRARY})

#link_libraries(/opt/xt-libsci/11.0.04.4/gnu/45/interlagos/lib/libsci_gnu.a)
