//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak. doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCHamiltonian/ReferencePoints.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_NEREFERENCE_POINTS_H
#define QMCPLUSPLUS_NEREFERENCE_POINTS_H

/** @file
 *  \todo When QMCHamiltonian/ReferencePoints is removed rename this class to ReferencePoints
 */

#include <Configuration.h>
#include "OhmmsData/OhmmsElementBase.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "OhmmsPETE/Tensor.h"
#include "ReferencePointsInput.h"

namespace qmcplusplus
{

/** This class creates, contains, and writes both user and machine readable referencepoints.
 *  they are derived from the lattice of the pset passed and the particle positions in the ref_psets
 *  an arbitrary number of additional points can be defined in the input that ReferencePointsInput
 *  presents as native input.
 *  It is a dependency of Estimators/NESpaceGrid and Estimatorss/EnergyDensityEstimator
 */
class NEReferencePoints
{
public:
  using Real   = QMCTraits::RealType;
  using Axes   = Tensor<Real, OHMMS_DIM>;
  using Point  = TinyVector<Real, OHMMS_DIM>;
  using Points = std::map<std::string, Point>;
  using Coord  = typename ReferencePointsInput::Coord;
  /** Usual constructor
   *  \param[in] rp_input   Input object for reference points which can contain and arbitrary set of points beyond
   *                        those take from the pset, and ref_psets
   *  \param[in] pset       pset that supplies the lattice information for reference points
   *  \param[in] ref_psets  pset reference vector the particle points in this/these psets are reference points
   */
  NEReferencePoints(const ReferencePointsInput& rp_input, const ParticleSet& pset, RefVector<ParticleSet>& ref_psets);
  NEReferencePoints(const NEReferencePoints& nerp) = default;

  /** writes a human readable representation of the reference points.
   *  \param[inout] os       ostream to write description to
   *  \param[in]    indent   spaces or other text to preface each line of output with.  needed to preserve
   *                         legacy output format.
   */
  void write_description(std::ostream& os, const std::string& indent) const;

  /** machine readable output
   *  \param[inout]  file  hdf5 file to write to.  Respects current state of file.
   */
  void write(hdf_archive& file) const;

  /** return const ref to map of reference points.
   *  labeling scheme unchanged from legacy
   */
  const Points& get_points() const { return points_; }

protected:
  Points points_;
private:
  void processParticleSets(const ParticleSet& P, RefVector<ParticleSet>& Pref);
  Axes axes;
  ReferencePointsInput input_;
};

std::ostream& operator<<(std::ostream& out, const NEReferencePoints& rhs);

}
#endif
