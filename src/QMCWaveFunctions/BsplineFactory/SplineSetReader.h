//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
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
#if defined(QMC_COMPLEX)
#include "SplineC2C.h"
#include "SplineC2COMPTarget.h"
#else
#include "SplineR2R.h"
#include "SplineC2R.h"
#include "SplineC2ROMPTarget.h"
#endif

namespace qmcplusplus
{
/** General SplineSetReader to handle any unitcell
 */
template<typename SA>
class SplineSetReader : public BsplineReader
{
public:
  using DataType   = typename SA::DataType;
  using SplineType = typename SA::SplineType;

  SplineSetReader(EinsplineSetBuilder* e);

  std::unique_ptr<SPOSet> create_spline_set(const std::string& my_name,
                                            int spin,
                                            const BandInfoGroup& bandgroup) override;

  /** create data space in the spline object and try open spline dump files.
   * @param bandgroup band info
   * @param bspline the spline object being worked on
   * @return true if dumpfile pass class name and data type size check
   */
  bool createSplineDataSpaceLookforDumpFile(const BandInfoGroup& bandgroup, SA& bspline) const;

  /** read planewave coefficients from h5 file
   * @param s data set full path in h5
   * @param h5f hdf5 file handle
   * @param cG vector to store coefficients
   */
  void readOneOrbitalCoefs(const std::string& s, hdf_archive& h5f, Vector<std::complex<double>>& cG) const;

  /** transforming planewave orbitals to 3D B-spline orbitals in real space.
   * @param spin orbital dataset spin index
   * @param bandgroup band info
   * @param bspline the spline object being worked on
   */
  void initialize_spline_pio_gather(const int spin, const BandInfoGroup& bandgroup, SA& bspline) const;
};

#if defined(QMC_COMPLEX)
extern template class SplineSetReader<SplineC2C<float>>;
extern template class SplineSetReader<SplineC2COMPTarget<float>>;
#if !defined(QMC_MIXED_PRECISION)
extern template class SplineSetReader<SplineC2C<double>>;
extern template class SplineSetReader<SplineC2COMPTarget<double>>;
#endif
#else
extern template class SplineSetReader<SplineC2R<float>>;
extern template class SplineSetReader<SplineR2R<float>>;
extern template class SplineSetReader<SplineC2ROMPTarget<float>>;
#if !defined(QMC_MIXED_PRECISION)
extern template class SplineSetReader<SplineC2R<double>>;
extern template class SplineSetReader<SplineR2R<double>>;
extern template class SplineSetReader<SplineC2ROMPTarget<double>>;
#endif
#endif

} // namespace qmcplusplus
#endif
