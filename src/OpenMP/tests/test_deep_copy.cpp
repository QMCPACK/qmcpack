//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <memory>
#include <vector>
#include <iostream>
#include "OpenMP/OMPallocator.hpp"

namespace qmcplusplus
{

struct container
{
  double* data;
  int size;
};

TEST_CASE("OMPdeepcopy", "[OMP]")
{
  const int MAX = 100;
  auto* foo = new container;

  OMPallocator<double> myAlloc;

  foo->data = myAlloc.allocate(MAX);
  for(int i=0; i<MAX; i++)
    foo->data[i] = i;
  foo->size = MAX;

  auto* data_ptr = foo->data;
  PRAGMA_OMP("omp target update to(data_ptr[0:foo->size])")

  int check_size(0);
  double check_data1(0);
  PRAGMA_OMP("omp target teams num_teams(1) map(from: check_size, check_data1)")
  {
     check_size = foo->size;
     check_data1 = foo->data[1];
  }

  REQUIRE(check_size == MAX);
  REQUIRE(check_data1 == 1.0);

  myAlloc.deallocate(foo->data,MAX);
  delete foo;
}

}
