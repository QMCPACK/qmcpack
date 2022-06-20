//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**@file CUDAruntime.hpp
 *@brief handle CUDA/HIP runtime selection.
 */
#ifndef QMCPLUSPLUS_CUDA_RUNTIME_H
#define QMCPLUSPLUS_CUDA_RUNTIME_H

#include <cstddef>
#include "config.h"
#ifndef QMC_CUDA2HIP
#include <cuda_runtime.h>
#else
#include <hip/hip_runtime.h>
#include "Platforms/ROCm/cuda2hip.h"
#endif

#include "CUDAerror.h"

size_t getCUDAdeviceFreeMem();

#endif
