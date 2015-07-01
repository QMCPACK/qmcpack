SET(CMAKE_SYSTEM_PROCESSOR "Haswell")

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)
set(GNU_OPTS "-DADD_ -DINLINE_ALL=inline")
set(GNU_FLAGS " -fomit-frame-pointer -malign-double -fopenmp -g -O3 -Drestrict=__restrict__ -finline-limit=1000 -fstrict-aliasing -funroll-all-loops ")
set(XT_FLAGS " -msse -msse2 -msse3 -msse4.1 ")
set(CMAKE_CXX_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -ftemplate-depth-60 ${GNU_OPTS} -Wno-deprecated")
set(CMAKE_C_FLAGS "${XT_FLAGS} ${GNU_FLAGS} -std=c99")

# need for both c++ and c
set(QMC_CUDA 1)
SET(ENABLE_OPENMP 1)
SET(HAVE_MPI 1)
set(HAVE_CUDA 1)
#SET(QMC_BUILD_STATIC 1)
message(STATUS "QMC_CUDA: " ${QMC_CUDA})

SET(HAVE_SSE 1)
SET(HAVE_SSE2 1)
SET(HAVE_SSE3 1)
SET(HAVE_SSSE3 1)
SET(HAVE_SSE41 1)
SET(USE_PREFETCH 1)
SET(PREFETCH_AHEAD 12)

SET(CMAKE_CXX_LINK_SHARED_LIBRARY)
SET(CMAKE_CXX_LINK_MODULE_LIBRARY)
SET(CMAKE_C_LINK_SHARED_LIBRARY)
SET(CMAKE_C_LINK_MODULE_LIBRARY)

set(CMAKE_FIND_ROOT_PATH
     /home/projects/qmcpack/boost_1_45_0
     /soft/libraries/hdf5-1.8.13
)

include_directories(/soft/compilers/intel/mkl/include)
link_libraries(-L/soft/compilers/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm)

# bypass einspline search
#set(EINSPLINE_HOME /lustre/widow3/scratch/jnkim/einspline)
#set(HAVE_EINSPLINE 1)
#set(HAVE_EINSPLINE_EXT 0)

set(CUDA_NVCC_FLAGS "-arch=sm_35;-Drestrict=__restrict__;-DNO_CUDA_MAIN;-O3")
set(CUDA_CUDART_LIBRARY /soft/visualization/cuda-7.0.28/lib64/libcudart.so)
#set(CUDA_CUDA_LIBRARY /soft/visualization/cuda-7.0.28/lib64/stubs/libcuda.so)
set(CUDA_cublas_LIBRARY /soft/visualization/cuda-7.0.28/lib64/libcublas.so)
set(CUDA_TOOLKIT_ROOT_DIR /soft/visualization/cuda-7.0.28)
set(CUDA_TOOLKIT_INCLUDE /soft/visualization/cuda-7.0.28/include)
set(CUDA_LIBRARIES ${CUDA_CUDART_LIBRARY})
set(CUDA_make2cmake ${CMAKE_ROOT}/Modules/FindCUDA/make2cmake.cmake)
set(CUDA_parse_cubin ${CMAKE_ROOT}/Modules/FindCUDA/parse_cubin.cmake)
set(CUDA_run_nvcc ${CMAKE_ROOT}/Modules/FindCUDA/run_nvcc.cmake)
set(CUDA_PROPAGATE_HOST_FLAGS OFF) #do not propagate the host flags


message(STATUS "CUDA_LIBRARY" $ENV{CUDA_CUDART_LIBRARY})
