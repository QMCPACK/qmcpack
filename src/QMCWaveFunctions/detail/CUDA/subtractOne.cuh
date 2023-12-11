//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef CUDA_SUBTRACTONE_H
#define CUDA_SUBTRACTONE_H

#include "config.h"
#ifndef QMC_CUDA2HIP
#include <cuComplex.h>
#else
#include <hip/hip_complex.h>
#include "ROCm/cuda2hip.h"
#endif

template<typename T>
__host__ __device__ __inline__ T subtractOne(T x)
{
  return x+T(-1);
}

template<>
__host__ __device__ __inline__ cuComplex subtractOne<cuComplex>(cuComplex x)
{
  return make_cuComplex(cuCrealf(x)-1.0f, cuCimagf(x));
}

template<>
__host__ __device__ __inline__ cuDoubleComplex subtractOne<cuDoubleComplex>(cuDoubleComplex x)
{
  return make_cuDoubleComplex(cuCreal(x)-1.0, cuCimag(x));
}

#endif
