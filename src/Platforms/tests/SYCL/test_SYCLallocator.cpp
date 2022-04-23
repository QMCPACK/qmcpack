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
#include <iostream>
#include "SYCL/SYCLallocator.hpp"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{
TEST_CASE("SYCL_allocator", "[SYCL]")
{
  // SYCLAllocator
  sycl::queue m_queue = getSYCLDefaultDeviceDefaultQueue();
  Vector<double, SYCLAllocator<double>> vec(1024);
  Vector<double> vec_h(1024);

  sycl::event e;
  {
    double* V = vec.data();

    e = m_queue.parallel_for(sycl::range<1>{1024}, [=](sycl::id<1> item) { V[item] = item + 1; });
  }

  //copy to host
  m_queue.memcpy(vec_h.data(), vec.data(), 1024 * sizeof(double), {e}).wait();

  CHECK(vec_h[0] == 1);
  CHECK(vec_h[77] == 78);
}

TEST_CASE("SYCL_host_allocator", "[SYCL]")
{
  sycl::queue m_queue = getSYCLDefaultDeviceDefaultQueue();
  // SYCLHostAllocator
  Vector<double, SYCLHostAllocator<double>> vec(1024, 1);

  {
    double* V = vec.data();
    m_queue.parallel_for(sycl::range<1>{1024}, [=](sycl::id<1> item) { V[item] += item + 1; }).wait();
  }

  CHECK(vec[0] == 2);
  CHECK(vec[77] == 79);
}

/*
TEST_CASE("SYCL_shared_allocator", "[SYCL]")
{
  sycl::queue m_queue = getSYCLDefaultDeviceDefaultQueue();
  Vector<double, SYCLSharedAllocator<double>> vec(1024);

  std::cout << "Size " << vec.size() << std::endl;
  {
    double* V = vec.data();
    m_queue.parallel_for(sycl::range<1>{1024}, [=](sycl::id<1> item) { V[item] = item + 1; }).wait();
  }
  CHECK(vec[0] == 1);
  CHECK(vec[77] == 78);
}
*/

} // namespace qmcplusplus
