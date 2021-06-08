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


#ifndef QMCPLUSPLUS_DTMODES_H
#define QMCPLUSPLUS_DTMODES_H

#include <cstdint>

namespace qmcplusplus
{
enum class DTModes : uint_fast8_t
{
  ALL_OFF = 0x0,
  /** whether full table needs to be ready at anytime or not
   * Optimization can be implemented during forward PbyP move when the full table is not needed all the time.
   * DT consumers should know if full table is needed or not and request via addTable.
   */
  NEED_FULL_TABLE_ANYTIME   = 0x1,
  NO_NEED_TEMP_DATA_ON_HOST = 0x2
};

constexpr bool operator&(DTModes x, DTModes y)
{
  return static_cast<uint_fast8_t>(x) & static_cast<uint_fast8_t>(y) != 0x0;
}

constexpr DTModes operator|(DTModes x, DTModes y)
{
  return static_cast<DTModes>(static_cast<uint_fast8_t>(x) | static_cast<uint_fast8_t>(y));
}

inline constexpr DTModes operator~(DTModes x) { return static_cast<DTModes>(~static_cast<uint_fast8_t>(x)); }

inline DTModes& operator|=(DTModes& x, DTModes y)
{
  x = x | y;
  return x;
}
} // namespace qmcplusplus
#endif
