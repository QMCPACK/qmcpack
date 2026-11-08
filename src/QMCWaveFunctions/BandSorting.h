//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

/** @file BandSorting.h
 * @brief Shared in-place band ordering for spline and plane-wave builders.
 */
#ifndef QMCPLUSPLUS_BAND_SORTING_H
#define QMCPLUSPLUS_BAND_SORTING_H

#include <algorithm>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace qmcplusplus
{
/** Order orbital bands, with deterministic treatment of near-degenerate energies.
 *
 * Mode 0 preserves input order. Mode 1 orders by energy; mode 2 first orders by
 * band index. Within each energy group, twist index then band index break ties.
 * Groups contain energies less than 1e-6 above their minimum, so a chain of
 * pairwise near-degenerate energies cannot join otherwise separate groups.
 * In mode 2, an energy group cannot span different band indices. All comparisons
 * passed to std::sort are strict lexicographic orderings; equal ordering keys
 * do not imply stable ordering of other band metadata.
 *
 * @param bands band records to reorder in place
 * @param sort_mode 0 for insertion order, 1 for energy, or 2 for band index
 * @param energy pointer to the double eigenvalue member of a band
 * @param twist_index pointer to the primitive k-point index member
 * @param band_index pointer to the band index member
 * @pre Eigenvalues are not NaN.
 * @throws std::runtime_error if sort_mode is not 0, 1, or 2.
 */
template<class Band>
void sortOrbitalBands(std::vector<Band>& bands,
                      int sort_mode,
                      double Band::* energy,
                      int Band::* twist_index,
                      int Band::* band_index)
{
  if (sort_mode == 0)
    return;
  if (sort_mode != 1 && sort_mode != 2)
    throw std::runtime_error("Orbital sort must be 0, 1, or 2");

  std::sort(bands.begin(), bands.end(), [&](const Band& left, const Band& right) {
    if (sort_mode == 2 && left.*band_index != right.*band_index)
      return left.*band_index < right.*band_index;
    return std::tie(left.*energy, left.*twist_index, left.*band_index) <
        std::tie(right.*energy, right.*twist_index, right.*band_index);
  });

  // Anchor each near-degenerate group at its minimum, rather than chaining tolerances.
  for (auto first = bands.begin(); first != bands.end();)
  {
    const double minimum_energy = (*first).*energy;
    const int group_band_index  = (*first).*band_index;
    auto last                   = first + 1;
    while (last != bands.end() && (sort_mode == 1 || (*last).*band_index == group_band_index) &&
           (*last).*energy - minimum_energy < 1e-6)
      ++last;
    std::sort(first, last, [&](const Band& left, const Band& right) {
      return std::tie(left.*twist_index, left.*band_index, left.*energy) <
          std::tie(right.*twist_index, right.*band_index, right.*energy);
    });
    first = last;
  }
}
} // namespace qmcplusplus

#endif
