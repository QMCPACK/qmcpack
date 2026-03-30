//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@inte.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * The most general reader class for the following classes using the full single grid for the supercell
 * - SplineR2R
 * - SplineC2C
 * - SplineC2R
 * Each band is initialized with UBspline_3d_d and both real and imaginary parts are passed to the objects
 * which will convert the data type to their internal precision. 
 */
#ifndef QMCPLUSPLUS_SPLINESET_READER_H
#define QMCPLUSPLUS_SPLINESET_READER_H

#include "BsplineReader.h"

namespace qmcplusplus
{

// forward declaration
template<typename ST>
class MultiBsplineBase;

/** General SplineSetReader to handle any unitcell
 */
template<typename ST>
class SplineSetReader : public BsplineReader
{
  std::unique_ptr<SPOSet> create_spline_set(const std::string& my_name,
                                            int spin,
                                            const BandInfoGroup& bandgroup) override;

public:
  SplineSetReader(EinsplineSetBuilder* e, bool use_duplex_splines);

  /** transforming planewave orbitals to 3D B-spline orbitals in real space.
   * @param spin orbital dataset spin index
   * @param bandgroup band info
   * @param bspline the spline object being worked on
   */
  void initialize_spline_pio_gather(const int spin,
                                    const BandInfoGroup& bandgroup,
                                    const TinyVector<int, 3>& half_g,
                                    const aligned_vector<int>& band_index_map,
                                    MultiBsplineBase<ST>& multi_splines) const;
};

extern template class SplineSetReader<float>;
#if !defined(QMC_MIXED_PRECISION)
extern template class SplineSetReader<double>;
#endif

} // namespace qmcplusplus
#endif
