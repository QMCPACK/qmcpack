//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak. doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCHamiltonian/ReferencePoints.h
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NEREFERENCE_POINTS_H
#define QMCPLUSPLUS_NEREFERENCE_POINTS_H

#include <Configuration.h>
#include "OhmmsData/OhmmsElementBase.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "OhmmsPETE/Tensor.h"
#include "ReferencePointsInput.h"

namespace qmcplusplus
{
class NEReferencePoints
{
public:
  using Real   = QMCTraits::RealType;
  using Axes   = Tensor<Real, OHMMS_DIM>;
  using Point  = TinyVector<Real, OHMMS_DIM>;
  using Points = std::map<std::string, Point>;
  using Coord  = typename ReferencePointsInput::Coord;
  NEReferencePoints(const ReferencePointsInput& rp_input, ParticleSet& pset, RefVector<ParticleSet>& ref_psets);
  NEReferencePoints(const NEReferencePoints& nerp) = default;

  void processParticleSets(ParticleSet& P, RefVector<ParticleSet>& Pref);
  void write_description(std::ostream& os, const std::string& indent) const;
  void write(hdf_archive& file) const;
  const Points& get_points() const { return points_; }

protected:
  Points points_;

private:
  Axes axes;
  ReferencePointsInput input_;
};

std::ostream& operator<<(std::ostream& out, const NEReferencePoints& rhs);

namespace testing
{
class TestableNEReferencePoints : public NEReferencePoints
{
public:
  TestableNEReferencePoints(const NEReferencePoints& nerp) : NEReferencePoints(nerp) {}
  void write_testable_description(std::ostream& os) const;
};
} // namespace testing

std::ostream& operator<<(std::ostream& out, const testing::TestableNEReferencePoints& rhs);
} // namespace qmcplusplus


#endif
