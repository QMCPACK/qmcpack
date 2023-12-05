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

set(CUDA_RT_SRCS CUDAfill.cpp CUDAallocator.cpp CUDAruntime.cpp CUDADeviceManager.cpp)
set(CUDA_LA_SRCS cuBLAS_missing_functions.cu)

add_library(platform_cuda_runtime ${CUDA_RT_SRCS})
add_library(platform_cuda_LA ${CUDA_LA_SRCS})

target_link_libraries(platform_cuda_runtime PRIVATE platform_host_runtime)

if(NOT QMC_CUDA2HIP)
  target_link_libraries(platform_cuda_runtime PUBLIC CUDA::cudart)
  target_link_libraries(platform_cuda_LA PUBLIC CUDA::cublas CUDA::cusolver)
else()
  set_source_files_properties(${CUDA_LA_SRCS} PROPERTIES LANGUAGE HIP)
  target_link_libraries(platform_cuda_runtime PUBLIC platform_rocm_runtime)
  target_link_libraries(platform_cuda_LA PRIVATE qmc_external_thrust)
  target_link_libraries(platform_cuda_LA PUBLIC platform_rocm_LA)
endif()

target_link_libraries(platform_cuda_LA PUBLIC platform_cuda_runtime)
