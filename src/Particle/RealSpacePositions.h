//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file RealSpacePostions.h
 */
#ifndef QMCPLUSPLUS_REALSPACE_POSITIONS_H
#define QMCPLUSPLUS_REALSPACE_POSITIONS_H

#include "Particle/DynamicCoordinates.h"
#include "OhmmsSoA/VectorSoaContainer.h"

namespace qmcplusplus
{
/** Introduced to handle virtual moves and ratio computations, e.g. for non-local PP evaluations.
   */
class RealSpacePositions : public DynamicCoordinates
{
public:
  using ParticlePos = PtclOnLatticeTraits::ParticlePos;
  using RealType    = QMCTraits::RealType;
  using PosType     = QMCTraits::PosType;

  RealSpacePositions() : DynamicCoordinates(DynamicCoordinateKind::DC_POS) {}

  std::unique_ptr<DynamicCoordinates> makeClone() override { return std::make_unique<RealSpacePositions>(*this); }

  void resize(size_t n) override { RSoA.resize(n); }
  size_t size() const override { return RSoA.size(); }

  void setAllParticlePos(const ParticlePos& R) override
  {
    resize(R.size());
    RSoA.copyIn(R);
  }
  void setOneParticlePos(const PosType& pos, size_t iat) override { RSoA(iat) = pos; }

  void mw_acceptParticlePos(const RefVectorWithLeader<DynamicCoordinates>& coords_list,
                            size_t iat,
                            const std::vector<PosType>& new_positions,
                            const std::vector<bool>& isAccepted) const override
  {
    assert(this == &coords_list.getLeader());
    for (size_t iw = 0; iw < isAccepted.size(); iw++)
      if (isAccepted[iw])
        coords_list[iw].setOneParticlePos(new_positions[iw], iat);
  }

  const PosVectorSoa& getAllParticlePos() const override { return RSoA; }
  PosType getOneParticlePos(size_t iat) const override { return RSoA[iat]; }

private:
  ///particle positions in SoA layout
  PosVectorSoa RSoA;
};
} // namespace qmcplusplus
#endif
