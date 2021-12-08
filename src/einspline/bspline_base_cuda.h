//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (C) 2007-2010 Kenneth P. Esler, Jr.
// Modifications Copyright (C) 2021 Advanced Micro Devices, Inc. All rights reserved.
//
// File developed by: Kenneth P. Esler, Jr.
//                    Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//
// File created by: Kenneth P. Esler, Jr.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef BSPLINE_BASE_CUDA_H
#define BSPLINE_BASE_CUDA_H

#include "config.h"
#ifndef QMC_CUDA2HIP
#include <cuda.h>
#else
#include <hip/hip_runtime.h>
#include "Platforms/ROCm/cuda2hip.h"
#endif

#if defined(CUDA_VERSION) && (CUDA_VERSION < 3000) /* 3.0 */
typedef struct
{
  double x,y,z;
} double3;

typedef struct
{
  double x,y,z,w;
} double4;
#endif

#endif
