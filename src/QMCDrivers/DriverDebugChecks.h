//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DRIVERDEBUGCHECKS_H
#define QMCPLUSPLUS_DRIVERDEBUGCHECKS_H

#include <cstdint>

namespace qmcplusplus
{
enum class DriverDebugChecks : uint_fast8_t
{
  ALL_OFF             = 0x0,
  CHECKGL_AFTER_LOAD  = 0x1,
  CHECKGL_AFTER_MOVES = 0x2,
  CHECKGL_AFTER_TMOVE = 0x3,
};

constexpr bool operator&(DriverDebugChecks x, DriverDebugChecks y)
{
  return (static_cast<uint_fast8_t>(x) & static_cast<uint_fast8_t>(y)) != 0x0;
}

constexpr DriverDebugChecks operator|(DriverDebugChecks x, DriverDebugChecks y)
{
  return static_cast<DriverDebugChecks>(static_cast<uint_fast8_t>(x) | static_cast<uint_fast8_t>(y));
}

constexpr DriverDebugChecks operator~(DriverDebugChecks x) { return static_cast<DriverDebugChecks>(~static_cast<uint_fast8_t>(x)); }

inline DriverDebugChecks& operator|=(DriverDebugChecks& x, DriverDebugChecks y)
{
  x = x | y;
  return x;
}
} // namespace qmcplusplus
#endif
