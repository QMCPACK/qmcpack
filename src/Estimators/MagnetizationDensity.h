//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MAGNETIZATION_DENSITY_H
#define QMCPLUSPLUS_MAGNETIZATION_DENSITY_H

#include <vector>
#include <functional>

#include "Estimators/OperatorEstBase.h"
#include "type_traits/complex_help.hpp"
#include "ParticleBase/RandomSeqGenerator.h"
#include <SpeciesSet.h>
#include <StdRandom.h>
#include "Estimators/MagnetizationDensityInput.h"

namespace qmcplusplus
{

class MagnetizationDensity : public OperatorEstBase
{
public:
  using Value         = QMCTraits::ValueType;
  using FullPrecValue = QMCTraits::FullPrecValueType;
  using Real          = RealAlias<Value>;
  using FullPrecReal  = RealAlias<FullPrecValue>;
  using Grad          = TinyVector<Value, OHMMS_DIM>;
  using Lattice       = PtclOnLatticeTraits::ParticleLayout;
  using Position      = QMCTraits::PosType;

  MagnetizationDensity(DataLocality dl):OperatorEstBase(dl){};

  void startBlock(int steps) override {};

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator& rng) override {};

  std::unique_ptr<OperatorEstBase> spawnCrowdClone() const override { return 0;};
private:
  MagnetizationDensityInput input_;
};
} // namespace qmcplusplus
#endif  /* QMCPLUSPLUS_MAGNETIZATION_DENSITY_H */
