//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "PlatformSelector.hpp"
#include <config.h>

namespace qmcplusplus
{
TEST_CASE("PlatformSelector", "[platform]")
{
  SECTION("CPU_OMPTARGET")
  {
#if defined(ENABLE_OFFLOAD)
    CHECK(CPUOMPTargetSelector::selectPlatform("yes") == PlatformKind::OMPTARGET);
    CHECK(CPUOMPTargetSelector::selectPlatform("") == PlatformKind::OMPTARGET);
#else
    CHECK(CPUOMPTargetSelector::selectPlatform("yes") == PlatformKind::CPU);
    CHECK(CPUOMPTargetSelector::selectPlatform("") == PlatformKind::CPU);
#endif
    CHECK(CPUOMPTargetSelector::selectPlatform("omptarget") == PlatformKind::OMPTARGET);
    CHECK(CPUOMPTargetSelector::selectPlatform("cpu") == PlatformKind::CPU);
    CHECK(CPUOMPTargetSelector::selectPlatform("no") == PlatformKind::CPU);
  }

  SECTION("CPU_OMPTARGET_CUDA")
  {
    using CPUOMPTargetCUDASelector = PlatformSelector<SelectorKind::CPU_OMPTARGET_CUDA>;
#if defined(ENABLE_CUDA)
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("yes") == PlatformKind::CUDA);
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("") == PlatformKind::CUDA);
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("cuda") == PlatformKind::CUDA);
#elif defined(ENABLE_OFFLOAD)
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("yes") == PlatformKind::OMPTARGET);
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("") == PlatformKind::OMPTARGET);
#else
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("yes") == PlatformKind::CPU);
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("") == PlatformKind::CPU);
#endif
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("omptarget") == PlatformKind::OMPTARGET);
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("cpu") == PlatformKind::CPU);
    CHECK(CPUOMPTargetCUDASelector::selectPlatform("no") == PlatformKind::CPU);
  }

  SECTION("CPU_OMPTARGET_SYCL")
  {
    using CPUOMPTargetSYCLSelector = PlatformSelector<SelectorKind::CPU_OMPTARGET_SYCL>;
#if defined(ENABLE_SYCL)
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("yes") == PlatformKind::SYCL);
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("") == PlatformKind::SYCL);
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("sycl") == PlatformKind::SYCL);
#elif defined(ENABLE_OFFLOAD)
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("yes") == PlatformKind::OMPTARGET);
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("") == PlatformKind::OMPTARGET);
#else
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("yes") == PlatformKind::CPU);
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("") == PlatformKind::CPU);
#endif
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("omptarget") == PlatformKind::OMPTARGET);
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("cpu") == PlatformKind::CPU);
    CHECK(CPUOMPTargetSYCLSelector::selectPlatform("no") == PlatformKind::CPU);
  }
}
} // namespace qmcplusplus
