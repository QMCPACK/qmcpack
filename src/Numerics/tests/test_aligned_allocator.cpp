//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <iostream>
#include "config.h"
#include "simd/allocator.hpp"

namespace qmcplusplus
{

TEST_CASE("Aligned allocator", "[numerics]")
{
  bool not_aligned;

  aligned_vector<float> a(311);
  std::cout << "address=" << a.data() << " require=" << (void *)(QMC_CLINE-1) << std::endl;
  not_aligned = (size_t)a.data() & (QMC_CLINE-1);
  REQUIRE( !not_aligned );
  a.resize(829);
  std::cout << "address=" << a.data() << " require=" << (void *)(QMC_CLINE-1) << std::endl;
  not_aligned = (size_t)a.data() & (QMC_CLINE-1);
  REQUIRE( !not_aligned );

  aligned_vector<double> b(311);
  std::cout << "address=" << b.data() << " require=" << (void *)(QMC_CLINE-1) << std::endl;
  not_aligned = (size_t)b.data() & (QMC_CLINE-1);
  REQUIRE( !not_aligned );
  b.resize(829);
  std::cout << "address=" << b.data() << " require=" << (void *)(QMC_CLINE-1) << std::endl;
  not_aligned = (size_t)b.data() & (QMC_CLINE-1);
  REQUIRE( !not_aligned );
}

}
