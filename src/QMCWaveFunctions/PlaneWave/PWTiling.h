//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

/** @file PWTiling.h
 * @brief Primitive-twist folding, time-reversal pairing, and lattice validation.
 *
 * These helpers operate on supplied metadata without a communicator. Validation
 * failures use std::runtime_error; the builder supplies communicator-level handling.
 */
#ifndef QMCPLUSPLUS_PW_TILING_H
#define QMCPLUSPLUS_PW_TILING_H

#include "Configuration.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/TinyVector.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace qmcplusplus::pw
{
/// Reduced-coordinate twist represented in full precision.
using Twist = TinyVector<double, OHMMS_DIM>;
/// Integer matrix mapping primitive-cell lattice vectors to the supercell.
using TileMatrix = Tensor<int, OHMMS_DIM>;
/// Full-precision real-space lattice matrix.
using Lattice = Tensor<double, OHMMS_DIM>;

/** Primitive twists grouped by the supercell twist to which they fold.
 *
 * Entries at the same index in super_twists and primitive_indices describe one
 * supercell twist and the primitive-cell twist indices that contribute to it.
 */
struct TwistGroups
{
  /// Canonical supercell twists, in first-encounter order from the input file.
  std::vector<Twist> super_twists;
  /// Primitive k-point indices for each group, preserving input-file order.
  std::vector<std::vector<int>> primitive_indices;
};

/** Representative primitive twist for constructing real-valued SPOs.
 * Self-inverse twists supply one orbital per band; paired twists can supply two.
 */
struct DistinctTwist
{
  /// Index of the retained member of a time-reversal pair.
  int twist_index;
  /// Whether the retained complex band supplies two real-valued SPOs.
  bool make_two_copies;
};

/** Map a reduced twist to the centered representative used by bspline input.
 *
 * Integer reciprocal-lattice translations are removed component by component,
 * producing the conventionally centered representative, including +0.5 at the
 * half-grid boundary.
 */
inline Twist canonicalizeTwist(Twist twist)
{
  for (int dimension = 0; dimension < OHMMS_DIM; ++dimension)
    twist[dimension] -= std::round(twist[dimension] - 1e-6);
  return twist;
}

/** Return whether two reduced twists differ only by a reciprocal-lattice vector.
 *
 * @param left first reduced-coordinate twist
 * @param right second reduced-coordinate twist
 * @param tolerance upper bound on the squared norm of the reduced difference
 */
inline bool equivalentTwists(const Twist& left, const Twist& right, double tolerance = 1e-6)
{
  Twist difference = left - right;
  for (int dimension = 0; dimension < OHMMS_DIM; ++dimension)
    difference[dimension] -= std::nearbyint(difference[dimension]);
  return dot(difference, difference) < tolerance;
}

/** Return whether two reduced twists form a time-reversal pair modulo integers.
 *
 * @param left first reduced-coordinate twist
 * @param right second reduced-coordinate twist
 * @param tolerance upper bound on the squared norm of the reduced sum
 */
inline bool timeReversedTwists(const Twist& left, const Twist& right, double tolerance = 1e-6)
{
  Twist sum = left + right;
  for (int dimension = 0; dimension < OHMMS_DIM; ++dimension)
    sum[dimension] -= std::nearbyint(sum[dimension]);
  return dot(sum, sum) < tolerance;
}

/** Select one representative from each time-reversal pair for a real build.
 *
 * Self-inverse twists are retained once. A non-self-inverse pair is represented
 * by its later entry in included_twists, marked to produce real and imaginary
 * SPO projections. Retained representatives preserve their input ordering.
 *
 * @param included_twists primitive k-point indices contributing to one supercell twist
 * @param primitive_twists reduced primitive k-points indexed by included_twists
 * @pre Included indices are valid and represent distinct primitive reciprocal classes.
 */
inline std::vector<DistinctTwist> findDistinctRealTwists(const std::vector<int>& included_twists,
                                                         const std::vector<Twist>& primitive_twists)
{
  std::vector<int> paired_twists;
  std::vector<int> distinct_twists;
  for (int i = 0; i < included_twists.size(); ++i)
  {
    bool distinct = true;
    for (int j = i + 1; j < included_twists.size(); ++j)
      if (timeReversedTwists(primitive_twists[included_twists[i]], primitive_twists[included_twists[j]]))
        distinct = false;
    if (distinct)
      distinct_twists.push_back(included_twists[i]);
    else
      paired_twists.push_back(included_twists[i]);
  }

  std::vector<DistinctTwist> result;
  result.reserve(distinct_twists.size());
  for (const int twist_index : distinct_twists)
  {
    const bool make_two_copies =
        std::find_if(paired_twists.begin(), paired_twists.end(), [&](const int paired_index) {
          return timeReversedTwists(primitive_twists[twist_index], primitive_twists[paired_index]);
        }) != paired_twists.end();
    result.push_back({twist_index, make_two_copies});
  }
  return result;
}

/** Convert an integer tile matrix to full-precision real storage. */
inline Lattice asRealMatrix(const TileMatrix& tile_matrix)
{
  Lattice result;
  for (int row = 0; row < OHMMS_DIM; ++row)
    for (int column = 0; column < OHMMS_DIM; ++column)
      result(row, column) = tile_matrix(row, column);
  return result;
}

/** Group primitive twists according to T k_primitive modulo integers.
 *
 * Each returned group contains primitive k-points that fold to the same
 * canonical supercell twist under tile_matrix T. Group indices follow the first
 * occurrence of each folded twist; they are not primitive k-point indices.
 * Grouping alone does not check completeness or duplicates: use validateTwistGroup().
 *
 * @param tile_matrix integer primitive-to-supercell transformation
 * @param primitive_twists reduced primitive k-points in orbital-file order
 */
inline TwistGroups groupPrimitiveTwists(const TileMatrix& tile_matrix, const std::vector<Twist>& primitive_twists)
{
  const Lattice tiling = asRealMatrix(tile_matrix);
  TwistGroups groups;
  for (int primitive_index = 0; primitive_index < primitive_twists.size(); ++primitive_index)
  {
    const Twist super_twist = canonicalizeTwist(dot(tiling, primitive_twists[primitive_index]));
    auto group              = std::find_if(groups.super_twists.begin(), groups.super_twists.end(),
                                           [&](const Twist& candidate) { return equivalentTwists(candidate, super_twist); });
    if (group == groups.super_twists.end())
    {
      groups.super_twists.push_back(super_twist);
      groups.primitive_indices.push_back({primitive_index});
    }
    else
      groups.primitive_indices[std::distance(groups.super_twists.begin(), group)].push_back(primitive_index);
  }
  return groups;
}

/** Find the group equivalent to a requested reduced supercell twist.
 *
 * @return zero-based folded-group index, matching groups.super_twists
 * @throws std::runtime_error if the requested twist is absent.
 */
inline int findTwistGroup(const TwistGroups& groups, const Twist& requested_twist)
{
  auto group = std::find_if(groups.super_twists.begin(), groups.super_twists.end(),
                            [&](const Twist& candidate) { return equivalentTwists(candidate, requested_twist); });
  if (group == groups.super_twists.end())
    throw std::runtime_error("Requested plane-wave supercell twist is not present in the HDF file");
  return std::distance(groups.super_twists.begin(), group);
}

/** Validate that a twist group is a complete, nonduplicated tiling set.
 *
 * A valid group contains abs(det(tile_matrix)) primitive twists with distinct
 * primitive reciprocal-lattice classes modulo integer translations. Duplicate
 * classes are compared after folding the reduced difference into supercell
 * coordinates, preserving distinct nearby classes in large tilings.
 *
 * @pre groups was formed from tile_matrix and primitive_twists with valid indices.
 * @throws std::runtime_error for a singular matrix, invalid group index,
 * incomplete group, or duplicate reciprocal coset.
 */
inline void validateTwistGroup(const TileMatrix& tile_matrix,
                               const std::vector<Twist>& primitive_twists,
                               const TwistGroups& groups,
                               int group_index)
{
  const int tiling_size = std::abs(det(tile_matrix));
  if (tiling_size == 0)
    throw std::runtime_error("Plane-wave tilematrix must be nonsingular");
  if (group_index < 0 || group_index >= groups.primitive_indices.size())
    throw std::runtime_error("Plane-wave supercell twist index is outside the available range");

  const std::vector<int>& primitive_indices = groups.primitive_indices[group_index];
  if (primitive_indices.size() != tiling_size)
    throw std::runtime_error("Plane-wave supercell twist does not contain abs(det(tilematrix)) primitive k-points");

  const Lattice tiling = asRealMatrix(tile_matrix);
  for (auto first = primitive_indices.begin(); first != primitive_indices.end(); ++first)
    for (auto second = first + 1; second != primitive_indices.end(); ++second)
    {
      // Remove primitive reciprocal translations before comparing. Measure the
      // residual in supercell coordinates, as in twist grouping, so nearby but
      // distinct primitive classes in large tilings are not treated as duplicates.
      const Twist difference = dot(tiling, canonicalizeTwist(primitive_twists[*first] - primitive_twists[*second]));
      if (dot(difference, difference) < 1e-6)
        throw std::runtime_error("Plane-wave supercell twist contains duplicate primitive k-point cosets");
    }
}

/** Construct the real-space supercell lattice T L_primitive. */
inline Lattice makeSuperLattice(const TileMatrix& tile_matrix, const Lattice& primitive_lattice)
{
  return dot(asRealMatrix(tile_matrix), primitive_lattice);
}

/** Verify that tilematrix and the primitive lattice reproduce the target lattice.
 *
 * @throws std::runtime_error when any nonzero matrix element differs beyond the
 * tolerance used by the bspline lattice check.
 */
inline void validateLattice(const TileMatrix& tile_matrix,
                            const Lattice& primitive_lattice,
                            const Lattice& target_lattice)
{
  const Lattice super_lattice = makeSuperLattice(tile_matrix, primitive_lattice);
  constexpr double tolerance  = 10 * std::numeric_limits<float>::epsilon();
  for (int row = 0; row < OHMMS_DIM; ++row)
    for (int column = 0; column < OHMMS_DIM; ++column)
    {
      const double scale = std::max(std::abs(super_lattice(row, column)), std::abs(target_lattice(row, column)));
      if (scale > tolerance && std::abs(super_lattice(row, column) - target_lattice(row, column)) / scale > tolerance)
        throw std::runtime_error("Plane-wave tilematrix and target lattice are inconsistent");
    }
}

/** Return whether a supercell twist can support real-valued boundary conditions.
 *
 * Each reduced component must be an integer or half-integer within the tolerance
 * used for the bspline-compatible check.
 */
inline bool isRealCompatible(const Twist& super_twist)
{
  constexpr double tolerance = 100 * 10 * std::numeric_limits<float>::epsilon();
  for (int dimension = 0; dimension < OHMMS_DIM; ++dimension)
    if (std::abs(2 * super_twist[dimension] - std::nearbyint(2 * super_twist[dimension])) > tolerance)
      return false;
  return true;
}
} // namespace qmcplusplus::pw

#endif
