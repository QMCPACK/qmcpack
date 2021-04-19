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
#include "Platforms/CUDA/CUDALinearAlgebraHandles.h"
#include "Platforms/PinnedAllocator.h"
#ifdef ENABLE_OFFLOAD
#include "OMPTarget/OMPallocator.hpp"
#elif ENABLE_CUDA
#include "DualAllocator.hpp"
#endif
#include "DeviceValue.hpp"
#include "test_kernel.hpp"

namespace qmcplusplus
{

#ifdef ENABLE_OFFLOAD
  template<typename T>
  using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;
#elif ENABLE_CUDA
  template<typename T>
  using OffloadPinnedAllocator = DualAllocator<T, CUDAAllocator<T>, PinnedAlignedAllocator<T>>;
#endif

TEST_CASE("DeviceValue simple access", "[Platforms]")
{
  auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();
  double value = 1.0;
  DeviceValue<double, OffloadPinnedAllocator<double>> test_value(value);
  bool result;
  cudaCheck(testing::checkValueCUDA(cuda_handles->hstream, test_value.getDevicePtr(), value, result));
  CHECK(result);
}
} // namespace qmcplusplus
