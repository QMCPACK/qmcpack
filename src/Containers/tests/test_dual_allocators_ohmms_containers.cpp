//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch.hpp>

#include <algorithm>
#include "config.h"
#include "PinnedAllocator.h"
#include "OMPTarget/OMPallocator.hpp"
#if defined(ENABLE_CUDA)
#include "DualAllocator.hpp"
#include "CUDA/CUDAallocator.hpp"
#elif defined(ENABLE_SYCL)
#include "DualAllocator.hpp"
#include "SYCL/SYCLallocator.hpp"
#endif
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "OhmmsSoA/VectorSoaContainer.h"

#include "makeRngSpdMatrix.hpp"
#include "Utilities/for_testing/checkMatrix.hpp"

namespace qmcplusplus
{
// This naming is pretty bad but consistent with the naming over much of the code
// its defined locally all over the place.
template<typename T>
using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;
#if defined(ENABLE_CUDA)
template<typename T>
using VendorDualPinnedAllocator = DualAllocator<T, CUDAAllocator<T>, PinnedAlignedAllocator<T>>;
#elif defined(ENABLE_SYCL)
template<typename T>
using VendorDualPinnedAllocator = DualAllocator<T, SYCLAllocator<T>, PinnedAlignedAllocator<T>>;
#endif

template<class OPA>
void testDualAllocator()
{
  constexpr int nd = 3;
  using Value      = typename OPA::value_type;
  VectorSoaContainer<Value, nd, OPA> vcsoa;

  int n = 12;
  Matrix<Value> mat_spd;
  mat_spd.resize(n, n);
  // Resize so each "dimension" can hold an n x n matrix
  vcsoa.resize(n * n);

  testing::MakeRngSpdMatrix<Value> makeRngSpdMatrix;

  for (int iw = 0; iw < nd; ++iw)
  {
    makeRngSpdMatrix(mat_spd);
    std::copy_n(mat_spd.data(), n * n, vcsoa.data(iw));
  }

  vcsoa.updateTo();

  Matrix<Value, OPA> matrix_view(vcsoa, vcsoa.data(0), n, n);
  Matrix<Value, OPA> matrix_view2(vcsoa, vcsoa.data(1), n, n);
  Matrix<Value, OPA> matrix_view3(vcsoa, vcsoa.data(2), n, n);


  // Without copying allocator this causes a segfault.
  auto device_ptr  = matrix_view.device_data();
  auto device_ptr2 = matrix_view2.device_data();
  auto device_ptr3 = matrix_view3.device_data();

  CHECK(device_ptr != nullptr);
  CHECK(device_ptr2 != nullptr);
  CHECK(device_ptr3 != nullptr);

  std::ptrdiff_t distance_host   = matrix_view.data() - vcsoa.data();
  std::ptrdiff_t distance_device = matrix_view.device_data() - vcsoa.device_data();
  CHECK(distance_host == distance_device);

  distance_host   = matrix_view2.data() - vcsoa.data();
  distance_device = matrix_view2.device_data() - vcsoa.device_data();
  CHECK(distance_host == distance_device);

  distance_host   = matrix_view3.data() - vcsoa.data();
  distance_device = matrix_view3.device_data() - vcsoa.device_data();
  CHECK(distance_host == distance_device);

  int ifrom = 2;
  int ito   = 0;

  vcsoa.copyDeviceDataByIndex(0, 2);
  vcsoa.updateFrom();

  auto check_matrix_result = checkMatrix(mat_spd, matrix_view);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  matrix_view2(0, 0) = 0.0;
  matrix_view2(1, 0) = 1, 0;

  matrix_view2.updateTo();
  vcsoa.copyDeviceDataByIndex(0, 1);
  matrix_view.updateFrom();
  CHECK(matrix_view(0, 0) == 0.0);
  CHECK(matrix_view(1, 0) == 1.0);

  Array<Value, 3, OPA> aa;
  aa.resize(2, 2, 3);
  CHECK(aa.size() == 12);
}

TEST_CASE("OhmmsMatrix_VectorSoaContainer_View", "[Integration][Allocators]")
{
  testDualAllocator<OffloadPinnedAllocator<double>>();
  testDualAllocator<OffloadPinnedAllocator<std::complex<double>>>();
#if defined(ENABLE_CUDA) || defined(ENABLE_SYCL)
  testDualAllocator<VendorDualPinnedAllocator<double>>();
  testDualAllocator<VendorDualPinnedAllocator<std::complex<double>>>();
#endif
}
} // namespace qmcplusplus
