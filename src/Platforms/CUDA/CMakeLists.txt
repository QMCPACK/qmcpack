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

set(CUDA_SRCS cuBLAS_missing_functions.cu CUDAfill.cpp CUDAallocator.cpp CUDAruntime.cpp)

if(NOT QMC_CUDA2HIP)
  cuda_add_library(platform_cuda ${CUDA_SRCS})
  cuda_add_cublas_to_target(platform_cuda)
  target_link_libraries(platform_cuda PRIVATE ${CUDA_LIBRARIES} ${CUDA_cusolver_LIBRARY})
else()
  hip_add_library(platform_cuda ${CUDA_SRCS})
  target_link_libraries(platform_cuda PUBLIC hip::host roc::hipblas roc::rocsolver)
  target_include_directories(platform_cuda PUBLIC ${PROJECT_SOURCE_DIR}/external_codes/thrust)
endif()
