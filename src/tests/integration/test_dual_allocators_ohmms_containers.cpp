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
#if defined(ENABLE_OFFLOAD)
#warning "OMP allocator included"
#include "OMPTarget/OMPallocator.hpp"
#endif
#include "DualAllocator.hpp"
#include "CUDA/CUDAallocator.hpp"

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsSoA/VectorSoaContainer.h"

#include "makeRngSpdMatrix.hpp"
#include "Utilities/for_testing/checkMatrix.hpp"

namespace qmcplusplus
{

#if defined(ENABLE_OFFLOAD)
template<typename T>
using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;
#elif defined(ENABLE_CUDA)
template<typename T>
using OffloadPinnedAllocator = DualAllocator<T, CUDAAllocator<T>, PinnedAlignedAllocator<T>>;
#endif


TEST_CASE("OhmmsMatrix_VectorSoaContainer_View", "[Integration][Allocators]")
{
  constexpr int nd      = 5;
  using Value = double;
  VectorSoaContainer<Value, nd, OffloadPinnedAllocator<Value>> vcsoa;

  int n = 12;
  Matrix<Value> mat_spd;
  mat_spd.resize(n, n);
  // Resize so each "dimension" can hold an n x n matrix
  vcsoa.resize(n * n);

  for (int iw = 0; iw < nd; ++iw)
  {
    testing::makeRngSpdMatrix(mat_spd);
    std::copy_n(vcsoa.data(iw), n * n, mat_spd.data());
  }

  vcsoa.updateTo();

  Matrix<Value, OffloadPinnedAllocator<Value>> matrix_view(vcsoa, vcsoa.data(0), n, n);
  Matrix<Value, OffloadPinnedAllocator<Value>> matrix_view2(vcsoa, vcsoa.data(1), n, n);

  // Without copying allocator this causes a segfault.
  auto device_ptr = matrix_view.device_data();
  auto device_ptr2 = matrix_view2.device_data();
  CHECK(device_ptr != nullptr);
  CHECK(device_ptr2 != nullptr);

  // more extensive testing that insures transfer is valid to and from would be a good idea
}
} // namespace qmcplusplus
