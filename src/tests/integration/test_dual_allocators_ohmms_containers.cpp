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
#include "OMPTarget/OMPallocator.hpp"
#endif
#if defined(ENABLE_CUDA)
#include "DualAllocator.hpp"
#include "CUDA/CUDAallocator.hpp"
#endif
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsSoA/VectorSoaContainer.h"

#include "makeRngSpdMatrix.hpp"
#include "Utilities/for_testing/checkMatrix.hpp"

namespace qmcplusplus
{

// This naming is pretty bad but consistent with the naming of OMPAllocator in the
// rest of the code.
#if defined(ENABLE_OFFLOAD)
template<typename T>
using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;
#endif
// To save confusion and allow testing both lets call this what it is.
#if defined(ENABLE_CUDA)
template<typename T>
using CUDAPinnedAllocator = DualAllocator<T, CUDAAllocator<T>, PinnedAlignedAllocator<T>>;
#endif

template<class OPA>
void testDualAllocator()
{
  constexpr int nd      = 3;
  using Value = typename OPA::value_type;
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
  auto device_ptr = matrix_view.device_data();
  auto device_ptr2 = matrix_view2.device_data();
  auto device_ptr3 = matrix_view3.device_data();

  CHECK(device_ptr != nullptr);
  CHECK(device_ptr2 != nullptr);
  CHECK(device_ptr3 != nullptr);

#if defined(ENABLE_CUDA)
  // not sure how to do this sort of thing device side with openmp target
  // cudaErrorCheck(cudaMemcpyAsync(device_ptr2, device_ptr, sizeof(Value) * n * n, cudaMemcpyDeviceToDevice),
  //                "cudaMemcpyAsync in test_dual_allocators");
  cudaErrorCheck(cudaMemcpyAsync(device_ptr, device_ptr3, sizeof(Value) * n * n, cudaMemcpyDeviceToDevice),
                 "cudaMemcpyAsync in test_dual_allocators");
  // cudaErrorCheck(cudaMemcpyAsync(device_ptr3, device_ptr2, sizeof(Value) * n * n, cudaMemcpyDeviceToDevice),
  //                "cudaMemcpyAsync in test_dual_allocators");
  cudaErrorCheck(cudaDeviceSynchronize(), "cudaStreamSynchronize failed!");

  vcsoa.updateFrom();

  auto check_matrix_result = checkMatrix(mat_spd, matrix_view);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
#endif
}

TEST_CASE("OhmmsMatrix_VectorSoaContainer_View", "[Integration][Allocators]")
{
#if defined(ENABLE_OFFLOAD)
  testDualAllocator<OffloadPinnedAllocator<double>>();
  //  testDualAllocator<OffloadPinnedAllocator<std::complex<double>>>();
#endif
#if defined(ENABLE_CUDA)
  testDualAllocator<CUDAPinnedAllocator<double>>();
  // testDualAllocator<CUDAPinnedAllocator<std::complex<double>>>();
#endif
}
} // namespace qmcplusplus
