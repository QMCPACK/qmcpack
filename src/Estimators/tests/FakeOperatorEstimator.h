//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_FAKEOPERATORESTIMATOR_H
#define QMCPLUSPLUS_FAKEOPERATORESTIMATOR_H

#include "OperatorEstBase.h"
#include "Configuration.h"
#include "type_traits/DataLocality.h"

namespace qmcplusplus
{
class FakeOperatorEstimator : public OperatorEstBase
{
public:
  using QMCT = QMCTraits;

  FakeOperatorEstimator(int num_ranks, DataLocality data_locality);

  FakeOperatorEstimator(const FakeOperatorEstimator& foe);

  ~FakeOperatorEstimator() override{};

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomBase<FullPrecRealType>& rng) override
  {}

  void registerOperatorEstimator(hdf_archive& file) override {}

  void startBlock(int nsteps) override {}

  std::unique_ptr<OperatorEstBase> spawnCrowdClone() const override
  {
    return std::make_unique<FakeOperatorEstimator>(*this);
  }

  void set_walker_weights(QMCT::RealType weight) { walkers_weight_ = weight; }
};

} // namespace qmcplusplus
#endif
