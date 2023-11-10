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
#ifndef QMCPLUSPLUS_REALSPACE_POSITIONST_H
#define QMCPLUSPLUS_REALSPACE_POSITIONST_H

#include "OhmmsSoA/VectorSoaContainer.h"
#include "Particle/DynamicCoordinatesT.h"

namespace qmcplusplus
{
/** Introduced to handle virtual moves and ratio computations, e.g. for
 * non-local PP evaluations.
 */
template<typename T>
class RealSpacePositionsT : public DynamicCoordinatesT<T>
{
public:
  using ParticlePos  = typename LatticeParticleTraits<T>::ParticlePos;
  using RealType     = typename DynamicCoordinatesT<T>::RealType;
  using PosType      = typename DynamicCoordinatesT<T>::PosType;
  using PosVectorSoa = typename DynamicCoordinatesT<T>::PosVectorSoa;

  RealSpacePositionsT() : DynamicCoordinatesT<T>(DynamicCoordinateKind::DC_POS) {}

  std::unique_ptr<DynamicCoordinatesT<T>> makeClone() override { return std::make_unique<RealSpacePositionsT>(*this); }

  void resize(size_t n) override { RSoA.resize(n); }
  size_t size() const override { return RSoA.size(); }

  void setAllParticlePos(const ParticlePos& R) override
  {
    resize(R.size());
    RSoA.copyIn(R);
  }
  void setOneParticlePos(const PosType& pos, size_t iat) override { RSoA(iat) = pos; }

  void mw_acceptParticlePos(const RefVectorWithLeader<DynamicCoordinatesT<T>>& coords_list,
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
  /// particle positions in SoA layout
  PosVectorSoa RSoA;
};
} // namespace qmcplusplus
#endif
