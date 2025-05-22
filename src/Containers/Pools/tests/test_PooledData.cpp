/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by:  Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Pools/PooledData.h"
#include "Utilities/Timer.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <random>
#include <complex>
#include <OhmmsPETE/OhmmsVector.h>


//#define CHECK_ALLOCATION_PERF

using std::complex;
using std::string;
using std::vector;
namespace qmcplusplus
{
TEST_CASE("PooledData::std::complex_into_real", "[utilities]")
{
  PooledData<double> pooled_data;
  Vector<std::complex<double>> vec_z{{1.0, 1.0}, {0.0, 1.0}, {2.0, -1.0}, {-3.0, 0.0}, {-1.0, -2.0}};
  pooled_data.add(vec_z.first_address(), vec_z.last_address());
  Vector<std::complex<double>> vec_z_get(5);
  pooled_data.rewind();
  pooled_data.get(vec_z_get.first_address(), vec_z_get.last_address());
  CHECK(vec_z == vec_z_get);
}

TEST_CASE("PooledData::reduced_prec_std::complex_into_real", "[utilities]")
{
  PooledData<double> pooled_data;
  Vector<std::complex<float>> vec_z{{1.0, 1.0}, {0.0, 1.0}, {2.0, -1.0}, {-3.0, 0.0}, {-1.0, -2.0}};
  pooled_data.add(vec_z.first_address(), vec_z.last_address());
  Vector<std::complex<float>> vec_z_get(5);
  pooled_data.rewind();
  pooled_data.get(vec_z_get.first_address(), vec_z_get.last_address());
  CHECK(vec_z == vec_z_get);
}


} // namespace qmcplusplus
