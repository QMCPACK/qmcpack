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
} // namespace qmcplusplus
