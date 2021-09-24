#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(CUDA_RT_SRCS CUDAfill.cpp CUDAallocator.cpp CUDAruntime.cpp)
set(CUDA_LA_SRCS cuBLAS_missing_functions.cu)

if(NOT QMC_CUDA2HIP)
  add_library(platform_cuda_runtime ${CUDA_RT_SRCS})
  add_library(platform_cuda_LA ${CUDA_LA_SRCS})
  target_link_libraries(platform_cuda_LA PRIVATE CUDA::cublas)
  target_include_directories(platform_cuda_runtime PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  target_link_libraries(platform_cuda_runtime PRIVATE CUDA::cudart)
else()
  hip_add_library(platform_cuda_runtime ${CUDA_RT_SRCS})
  hip_add_library(platform_cuda_LA ${CUDA_LA_SRCS})
  target_link_libraries(platform_cuda_runtime PUBLIC platform_rocm_runtime)
  target_link_libraries(platform_cuda_LA PUBLIC platform_rocm_LA)
endif()

target_link_libraries(platform_cuda_LA PRIVATE ${CUDA_cusolver_LIBRARY} platform_cuda_runtime)
