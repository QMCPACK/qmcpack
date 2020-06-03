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
#include "omp.h"
#include "OMPallocator.hpp"

namespace qmcplusplus
{
template<typename T>
struct maptest
{
  const static size_t size = 6;
  T data[size];

  PRAGMA_OFFLOAD("omp declare target")
  inline void set_value(const T* restrict array, size_t size, T* restrict data)
  {
    PRAGMA_OFFLOAD("omp for")
    for (int i = 0; i < size; i++)
      data[i] = array[i];
  }
  PRAGMA_OFFLOAD("omp end declare target")

  void run()
  {
    PRAGMA_OFFLOAD("omp target enter data map(alloc:data[0:6])")
    T newdata[size] = {0, 1, 2, 3, 4, 5};
    PRAGMA_OFFLOAD("omp target map(to:newdata[0:6]) map(from:data[0:6])")
    {
      PRAGMA_OFFLOAD("omp parallel")
      set_value(newdata, size, data);
    }
    PRAGMA_OFFLOAD("omp target exit data map(delete:data[0:6])")
    std::cout << "data[5] = " << data[5] << std::endl;
  }
};

TEST_CASE("OMPclass member", "[OMP]")
{
  maptest<double> tester;
  tester.run();
}

} // namespace qmcplusplus
