//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "einspline/bspline_base.h"
#include "einspline/bspline_create.h"
#include "einspline/bspline_eval_d.h"
#include "einspline/multi_bspline_create.h"
#include "einspline/multi_bspline_eval_d.h"
#include "einspline/multi_bspline_eval_s.h"

#include "QMCWaveFunctions/BsplineFactory/BsplineDeviceCUDA.h"

#include <iostream>

namespace qmcplusplus
{

TEST_CASE("BsplineDeviceCUDA_Instantiation", "[wavefunction]")
{
  BsplineDeviceCUDA<double, 3> BDC;
}

}

int bspline_force_cuda_link();

TEST_CASE("bspline_force_cuda_link", "[wavefunction]")
{
  int a = bspline_force_cuda_link();
}

