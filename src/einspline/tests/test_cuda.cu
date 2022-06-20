#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////


#include "einspline/bspline_base.h"
#include "einspline/bspline_create.h"
#include "einspline/bspline_create_cuda.h"
#include "einspline/bspline_eval_d.h"
#include "einspline/bspline_eval_cuda.h"
#include "einspline/multi_bspline_structs.h"
#include "einspline/multi_bspline_create.h"
#include "einspline/multi_bspline_create_cuda.h"
#include "einspline/multi_bspline_eval_cuda.h"
#include "einspline/multi_bspline_eval_s.h"
#include "Platforms/CUDA_legacy/cuda_error.h"
#include <cassert>


void test_multi(multi_UBspline_1d_s *cpuSpline, float *pos, float *vals_cuda)
{
  multi_UBspline_1d_s_cuda *gpuSpline = create_multi_UBspline_1d_s_cuda(cpuSpline);

  float *vals[1];
  float *pos_d;
  float **vals_d;

  cudaCheck(cudaMalloc((void**)&(pos_d),     sizeof(float)));
  cudaCheck(cudaMemcpy(pos_d,  pos,  sizeof(float), cudaMemcpyHostToDevice));


  cudaCheck(cudaMalloc((void**)&(vals_d),  sizeof(float*)));
  float *valBlock_d;
  cudaCheck(cudaMalloc((void**)&(valBlock_d),    sizeof(float)));
  vals[0] = valBlock_d;

  cudaCheck(cudaMemcpy(vals_d,  vals,  sizeof(float*), cudaMemcpyHostToDevice));

  eval_multi_multi_UBspline_1d_s_cuda(gpuSpline, pos_d, vals_d, 1);

  cudaCheck(cudaMemcpy(vals_cuda,  valBlock_d,  sizeof(float), cudaMemcpyDeviceToHost));
}
