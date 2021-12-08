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

set(CUDA_LEGACY_SRCS cuda_inverse.cu gpu_vector.cpp gpu_misc.cpp)
add_library(platform_cuda_legacy ${CUDA_LEGACY_SRCS})
if(NOT QMC_CUDA2HIP)
  target_link_libraries(platform_cuda_legacy PUBLIC CUDA::cublas CUDA::cudart)
else()
  set_source_files_properties(${CUDA_LEGACY_SRCS} PROPERTIES LANGUAGE HIP)
  target_link_libraries(platform_cuda_legacy PUBLIC platform_rocm_runtime platform_rocm_LA)
endif(NOT QMC_CUDA2HIP)
