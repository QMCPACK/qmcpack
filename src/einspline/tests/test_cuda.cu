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


#include "catch.hpp"

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

int force_cuda_link()
{
    return 1;
}

void test_multi(multi_UBspline_1d_s *cpuSpline, float *pos, float *vals_cuda)
{
  multi_UBspline_1d_s_cuda *gpuSpline = create_multi_UBspline_1d_s_cuda(cpuSpline);

  float *vals[1];
  float *pos_d;
  float **vals_d;

  cudaMalloc((void**)&(pos_d),     sizeof(float));
  REQUIRE(pos_d != NULL);
  cudaMemcpy(pos_d,  pos,  sizeof(float), cudaMemcpyHostToDevice);


  cudaMalloc((void**)&(vals_d),  sizeof(float*));
  float *valBlock_d;
  cudaMalloc((void**)&(valBlock_d),    sizeof(float));
  vals[0] = valBlock_d;

  cudaMemcpy(vals_d,  vals,  sizeof(float*), cudaMemcpyHostToDevice);

  eval_multi_multi_UBspline_1d_s_cuda(gpuSpline, pos_d, vals_d, 1);

  cudaMemcpy(vals_cuda,  valBlock_d,  sizeof(float), cudaMemcpyDeviceToHost);
}

TEST_CASE("multi_cuda_wrapper", "[einspline]")
{
  Ugrid x_grid;
  // GPU versions require the grid to start at zero
  x_grid.start = 0.0;
  x_grid.end = 10.0;
  x_grid.num = 2;

  float data[2];
  data[0] = 2.0;
  data[1] = 3.0;

  BCtype_s xBC;
  xBC.lCode = NATURAL;
  xBC.rCode = NATURAL;
  multi_UBspline_1d_s* s = create_multi_UBspline_1d_s(x_grid, xBC, 1);
  REQUIRE(s);
  set_multi_UBspline_1d_s(s, 0, data);

  float pos[1];
  pos[0] = 0.0;

  // Check the CPU value
  float cpu_val[1];
  eval_multi_UBspline_1d_s(s, pos[0], cpu_val);
  REQUIRE(cpu_val[0] == 2.0);

  // Check the GPU value
  float vals_output[1];
  test_multi(s, pos, vals_output);
  REQUIRE(vals_output[0] == 2.0);
}
