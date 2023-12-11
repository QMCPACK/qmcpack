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


#include "catch.hpp"

#include <iostream>
#include "config.h"
#include "MemoryUsage.h"

namespace qmcplusplus
{
TEST_CASE("OMP runtime memory", "[OMP]")
{
  PRAGMA_OFFLOAD("omp target")
  {
    // intentional empty target to initialize offload runtime library.
  }

  print_mem("OMP runtime memory", std::cout);
}

} // namespace qmcplusplus
