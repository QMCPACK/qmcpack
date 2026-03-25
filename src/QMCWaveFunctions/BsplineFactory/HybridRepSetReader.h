//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * derived from BsplineReader
 */

#ifndef QMCPLUSPLUS_HYBRIDREP_READER_H
#define QMCPLUSPLUS_HYBRIDREP_READER_H

#include "BsplineReader.h"

namespace qmcplusplus
{
// forward declaration
template<typename ST>
class MultiBsplineBase;
template<typename ST>
class HybridRepCenterOrbitals;

/** General HybridRepSetReader to handle any unitcell
 */
template<typename ST>
class HybridRepSetReader : public BsplineReader
{
  using HybridBase = HybridRepCenterOrbitals<ST>;

  std::unique_ptr<SPOSet> create_spline_set(const std::string& my_name,
                                            int spin,
                                            const BandInfoGroup& bandgroup) override;

public:
  HybridRepSetReader(EinsplineSetBuilder* e, bool use_duplex_splines);

  /** initialize basic parameters of atomic orbitals */
  void initialize_hybridrep_atomic_centers(HybridBase& multi_atomic_splines) const;

  /** initialize construct atomic orbital radial functions from plane waves */
  inline void create_atomic_centers_Gspace(const Vector<std::complex<double>>& cG,
                                           Communicate& band_group_comm,
                                           const int iorb,
                                           const std::complex<double>& rotate_phase,
                                           const TinyVector<int, 3>& half_g,
                                           HybridBase& multi_atomic_splines) const;

  /** transforming planewave orbitals to 3D B-spline orbitals and 1D B-spline radial orbitals in real space.
   * @param spin orbital dataset spin index
   * @param bandgroup band info
   * @param bspline the spline object being worked on
   */
  void initialize_hybrid_pio_gather(const int spin,
                                    const BandInfoGroup& bandgroup,
                                    const TinyVector<int, 3>& half_g,
                                    const aligned_vector<int>& band_index_map,
                                    HybridBase& multi_atomic_splines,
                                    MultiBsplineBase<ST>& multi_splines) const;
};

extern template class HybridRepSetReader<float>;
#if !defined(QMC_MIXED_PRECISION)
extern template class HybridRepSetReader<double>;
#endif
} // namespace qmcplusplus
#endif
