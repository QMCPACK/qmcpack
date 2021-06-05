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


SET(CUDA_LEGACY_SRCS
  cuda_inverse.cu
  gpu_vector.cpp
  gpu_misc.cpp
)

IF(NOT QMC_CUDA2HIP)
  CUDA_ADD_LIBRARY(platform_cuda_legacy ${CUDA_LEGACY_SRCS})
  CUDA_ADD_CUBLAS_TO_TARGET(platform_cuda_legacy)
ELSE()
  HIP_ADD_LIBRARY(platform_cuda_legacy ${CUDA_LEGACY_SRCS})
ENDIF(NOT QMC_CUDA2HIP)
TARGET_LINK_LIBRARIES(platform_cuda_legacy PRIVATE ${CUDA_LIBRARIES})
