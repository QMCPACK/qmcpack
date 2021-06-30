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
#include <omp.h>
#include "OMPTarget/OMPallocator.hpp"

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
  auto* foo     = new container;
  foo->size     = MAX;

  OMPallocator<double> myAlloc;

  foo->data = myAlloc.allocate(MAX);
  for (int i = 0; i < MAX; i++)
    foo->data[i] = i;

  auto* data_ptr = foo->data;

  PRAGMA_OFFLOAD("omp target enter data map(alloc:foo[0:1])")
  PRAGMA_OFFLOAD("omp target map(always, to: foo[0:1], data_ptr[0:foo->size])") { foo->data = data_ptr; }

  int check_size(0);
  double check_data1(0);
  void* check_address1(nullptr);
  void* check_address2(nullptr);
  void* check_address3(nullptr);
  // clang-format off
  PRAGMA_OFFLOAD("omp target teams num_teams(1) \
                  map(from: check_size, check_data1, check_address1, check_address2, check_address3)")
  // clang-format on
  {
    check_size     = foo->size;
    check_data1    = foo->data[1];
    check_address1 = data_ptr;
    check_address2 = foo->data;
    check_address3 = foo;
  }

  std::cout << "foo->data value on the host " << foo->data << std::endl;
  std::cout << "foo->data value on the device " << check_address2 << std::endl;
  std::cout << "foo->data mapped address on the device " << check_address1 << std::endl;
  std::cout << "foo value on the host " << foo << std::endl;
  std::cout << "foo mapped address on the device " << check_address3 << std::endl;

  REQUIRE(check_data1 == 1.0);
  REQUIRE(check_size == MAX);

  PRAGMA_OFFLOAD("omp target teams num_teams(1) map(always,from:data_ptr[0:foo->size])") { data_ptr[1] = 2; }

  REQUIRE(data_ptr[1] == 2.0);

  myAlloc.deallocate(foo->data, MAX);
  PRAGMA_OFFLOAD("omp target exit data map(delete:foo[0:1])")
  delete foo;
}

} // namespace qmcplusplus
