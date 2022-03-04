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
  CHECK(OMPTargetSelector::convertToPlatform("yes") == PlatformKind::OMPTARGET);
  CHECK(OMPTargetSelector::convertToPlatform("") == PlatformKind::OMPTARGET);
#else
  CHECK(OMPTargetSelector::convertToPlatform("yes") == PlatformKind::CPU);
  CHECK(OMPTargetSelector::convertToPlatform("") == PlatformKind::CPU);
#endif
  CHECK(OMPTargetSelector::convertToPlatform("omptarget") == PlatformKind::OMPTARGET);
  CHECK(OMPTargetSelector::convertToPlatform("cpu") == PlatformKind::CPU);
  CHECK(OMPTargetSelector::convertToPlatform("no") == PlatformKind::CPU);
}
}
