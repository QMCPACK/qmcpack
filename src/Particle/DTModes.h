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
  /** whether full table needs to be ready at anytime or not during PbyP
   * Optimization can be implemented during forward PbyP move when the full table is not needed all the time.
   * DT consumers should know if full table is needed or not and request via addTable.
   */
  NEED_FULL_TABLE_ANYTIME = 0x1,
  /** For distance tables of virtual particle (VP) sets constructed based on this table, whether full table is needed on host
   * The corresponding DT of VP need to set MW_EVALUATE_RESULT_NO_TRANSFER_TO_HOST accordingly.
   */
  NEED_VP_FULL_TABLE_ON_HOST = 0x2,
  /** whether temporary data set on the host is updated or not when a move is proposed.
   * Considering transferring data from accelerator to host is relatively expensive,
   * only request this when data on host is needed for unoptimized code path.
   * This flag affects three subroutines mw_move, mw_updatePartial, mw_finalizePbyP in DistanceTable.
   */
  NEED_TEMP_DATA_ON_HOST = 0x4,
  /** skip data transfer back to host after mw_evalaute full distance table.
   * this optimization can be used for distance table consumed directly on the device without copying back to the host.
   */
  MW_EVALUATE_RESULT_NO_TRANSFER_TO_HOST = 0x8,
  /** whether full table needs to be ready at anytime or not after donePbyP
   * Optimization can be implemented during forward PbyP move when the full table is not needed all the time.
   * DT consumers should know if full table is needed or not and request via addTable.
   */
  NEED_FULL_TABLE_ON_HOST_AFTER_DONEPBYP = 0x16,
};

constexpr bool operator&(DTModes x, DTModes y)
{
  return (static_cast<uint_fast8_t>(x) & static_cast<uint_fast8_t>(y)) != 0x0;
}

constexpr DTModes operator|(DTModes x, DTModes y)
{
  return static_cast<DTModes>(static_cast<uint_fast8_t>(x) | static_cast<uint_fast8_t>(y));
}

constexpr DTModes operator~(DTModes x) { return static_cast<DTModes>(~static_cast<uint_fast8_t>(x)); }

inline DTModes& operator|=(DTModes& x, DTModes y)
{
  x = x | y;
  return x;
}
} // namespace qmcplusplus
#endif
