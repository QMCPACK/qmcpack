//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

/** @file PWBandInfo.h
 * @brief Primitive-band descriptors and ordering for plane-wave SPO construction.
 */
#ifndef QMCPLUSPLUS_PW_BAND_INFO_H
#define QMCPLUSPLUS_PW_BAND_INFO_H

#include "QMCWaveFunctions/BandSorting.h"

namespace qmcplusplus::pw
{
/** Descriptor for one primitive-cell band used to construct a plane-wave SPO.
 *
 * Indices refer to the orbital file, within the spin channel selected by the
 * builder. In real builds, make_two_copies marks a retained time-reversal
 * representative that supplies real and imaginary projections;
 * use_imaginary_part distinguishes the imaginary entry after expansion.
 * Complex builds use one entry per complex band with both flags false.
 */
struct BandInfo
{
  /// Eigenvalue used by energy ordering.
  double energy;
  /// Zero-based primitive-cell twist index in the orbital HDF file.
  int twist_index;
  /// Zero-based band index within the primitive-cell k-point and spin channel.
  int band_index;
  /// Select the imaginary rather than real projection of a complex band.
  bool use_imaginary_part{false};
  /// Expand this complex band into independent real and imaginary SPOs.
  bool make_two_copies{false};
};

/** Apply the requested bspline-compatible ordering to a band list.
 *
 * sort_mode=0 preserves insertion order, sort_mode=1 orders by energy, and
 * sort_mode=2 orders primarily by band index. Near-degenerate groups use the
 * shared sortOrbitalBands rules; projection flags do not affect ordering.
 * This routine neither pairs twists nor expands two-copy bands.
 *
 * @param bands primitive-band descriptors to reorder in place
 * @param sort_mode ordering mode, matching the spline builder
 * @pre Eigenvalues are not NaN.
 * @throws std::runtime_error if sort_mode is not 0, 1, or 2.
 */
inline void sortBands(std::vector<BandInfo>& bands, int sort_mode)
{
  sortOrbitalBands(bands, sort_mode, &BandInfo::energy, &BandInfo::twist_index, &BandInfo::band_index);
}
} // namespace qmcplusplus::pw

#endif
