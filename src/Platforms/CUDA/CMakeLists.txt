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


SET(CUDA_SRCS
  cuBLAS_missing_functions.cu
  CUDAfill.cpp
  CUDAallocator.cpp
  CUDAruntime.cpp
)

CUDA_ADD_LIBRARY(platform_cuda ${CUDA_SRCS})
CUDA_ADD_CUBLAS_TO_TARGET(platform_cuda)
TARGET_LINK_LIBRARIES(platform_cuda PRIVATE ${CUDA_LIBRARIES} ${CUDA_cusolver_LIBRARY})
