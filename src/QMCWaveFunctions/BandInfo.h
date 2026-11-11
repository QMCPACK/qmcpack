//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file BandInfo.h
 * @brief Orbital-band metadata and band-group selection.
 */
#ifndef QMCPLUSPLUS_BANDINFO_H
#define QMCPLUSPLUS_BANDINFO_H
#include <Configuration.h>

namespace qmcplusplus
{
class SPOSetInfo;

/** Metadata for one primitive-cell band used by the spline builders.
 *
 * Ordering is supplied explicitly by sortOrbitalBands, rather than by a
 * comparison operator on this record. MakeTwoCopies describes expansion into
 * real orbitals, not a second entry in the primitive-band list.
 */
struct BandInfo
{
  /// Zero-based primitive k-point index in the orbital file.
  int TwistIndex;
  /// Zero-based band index within a primitive k-point and spin channel.
  int BandIndex;
  /// Band-group metadata.
  int BandGroup;
  /// Spin-channel metadata.
  int Spin;
  /// Single-particle eigenvalue used for band ordering.
  double Energy;
  /// Whether this complex band can supply two real SPOs from a time-reversal pair.
  bool MakeTwoCopies;
  /// Initialize an unselected band with BandIndex=-1 and a large energy sentinel.
  BandInfo() : TwistIndex(0), BandIndex(-1), BandGroup(0), Spin(0), Energy(1e9), MakeTwoCopies(false) {}
};

/** Consecutive entries selected from an already ordered primitive-band list.
 *
 * A band marked MakeTwoCopies can supply two SPOs, so the number of selected
 * band entries need not equal NumSPOs.
 */
struct BandInfoGroup
{
  /// Identifier assigned by the band-group consumer.
  int GroupID;
  /// Number of SPOs represented by the selected entries.
  int NumSPOs;
  /// Starting entry in the ordered primitive-band list.
  int FirstBand;
  /// Selected primitive-band entries, retaining the input order.
  std::vector<BandInfo> myBands;
  /// Name assigned to this group.
  std::string myName;
  /// Initialize an empty selection; the consumer assigns GroupID.
  BandInfoGroup();
  /// Return the number of selected primitive-band entries before SPO expansion.
  inline int getNumDistinctOrbitals() const { return myBands.size(); }
  /// Return the number of SPOs after accounting for two-copy bands.
  inline int getNumSPOs() const { return NumSPOs; }

  /** Replace the selection with consecutive bands supplying num_spos SPOs.
   * @param bigspace already ordered primitive-band list
   * @param first_orb starting entry in bigspace, not an expanded SPO index
   * @param num_spos number of SPOs to supply, including real-orbital copies
   * @pre first_orb is nonnegative and num_spos is positive.
   * @throws UniformCommunicateError if the starting entry is past the list or
   * the remaining bands cannot supply the requested SPO count.
   */
  void selectBands(const std::vector<BandInfo>& bigspace, int first_orb, int num_spos);
};

} // namespace qmcplusplus

#endif
