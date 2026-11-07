//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_PW_BAND_INFO_H
#define QMCPLUSPLUS_PW_BAND_INFO_H

#include "QMCWaveFunctions/BandSorting.h"

namespace qmcplusplus::pw
{
/** Descriptor for one primitive-cell band used to construct a plane-wave SPO.
 *
 * The projection flags describe how a complex band is expanded when QMCPACK is
 * built for real-valued wave functions.
 */
struct BandInfo
{
  /// Eigenvalue used by energy ordering.
  double energy;
  /// Primitive-cell twist index in the orbital HDF file.
  int twist_index;
  /// Band index within the primitive-cell twist and spin channel.
  int band_index;
  /// Select the imaginary rather than real projection of a complex band.
  bool use_imaginary_part{false};
  /// Expand this complex band into independent real and imaginary SPOs.
  bool make_two_copies{false};
};

/** Apply the requested bspline-compatible ordering to a band list.
 *
 * sort_mode=0 preserves insertion order, sort_mode=1 orders by energy, and
 * sort_mode=2 orders primarily by band index.
 *
 * @throws std::runtime_error if sort_mode is not 0, 1, or 2.
 */
inline void sortBands(std::vector<BandInfo>& bands, int sort_mode)
{
  sortOrbitalBands(bands, sort_mode, &BandInfo::energy, &BandInfo::twist_index, &BandInfo::band_index);
}
} // namespace qmcplusplus::pw

#endif
