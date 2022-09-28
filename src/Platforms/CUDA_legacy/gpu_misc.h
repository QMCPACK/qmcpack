//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
// Modifications Copyright (C) 2021 Advanced Micro Devices, Inc. All rights reserved.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef GPU_MISC_H
#define GPU_MISC_H

#include <cstdlib>
#include <cstdio>
#include <vector>

#include "config.h"
#include "cuda_error.h"
#ifndef QMC_CUDA2HIP
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#else
#include <hip/hip_runtime.h>
#include <hipblas.h>
#include "Platforms/ROCm/cuda2hip.h"
#endif

namespace gpu
{
extern cudaStream_t kernelStream;
extern cudaStream_t memoryStream;

extern cudaEvent_t syncEvent;

extern cudaEvent_t gradientSyncDiracEvent;
extern cudaEvent_t gradientSyncOneBodyEvent;
extern cudaEvent_t gradientSyncTwoBodyEvent;

extern cudaEvent_t ratioSyncDiracEvent;
extern cudaEvent_t ratioSyncOneBodyEvent;
extern cudaEvent_t ratioSyncTwoBodyEvent;

extern cublasHandle_t cublasHandle;

extern size_t MaxGPUSpineSizeMB;
extern int rank;
extern int relative_rank;
extern int device_group_size;
extern bool cudamps;
extern std::vector<int> device_group_numbers;
extern std::vector<int> device_rank_numbers;

void initCUDAStreams();
void initCUDAEvents();
void initCublas();

void finalizeCUDAStreams();
void finalizeCUDAEvents();
void finalizeCublas();

void synchronize();

void streamsSynchronize();

} // namespace gpu
#endif
