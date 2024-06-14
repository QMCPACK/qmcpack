//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "EstimatorBuffer.hpp"

namespace qmcplusplus
{

TEST_CASE("EstimatorBuffer::copyIn", "[estimators]")
{
  std::vector<int> cats(3);
  std::vector<int> dogs(2);
  std::vector<int> elephants(1);

  EstimatorBuffer<int> e_buff;
  auto size = e_buff.copyIn(1,cats,dogs,elephants);
  CHECK(size == 6);
  
}

}
