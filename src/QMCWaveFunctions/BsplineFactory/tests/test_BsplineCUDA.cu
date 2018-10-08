#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
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

int bspline_force_cuda_link()
{
    return 1;
}

TEST_CASE("bspline_maka", "[wavefunction]")
{
  Ugrid x_grid;
  // GPU versions require the grid to start at zero
  x_grid.start = 0.0;
  x_grid.end = 10.0;
  x_grid.num = 2;
}