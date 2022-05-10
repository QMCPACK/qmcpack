//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerBase.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "EstimatorManagerCrowd.h"

namespace qmcplusplus
{
EstimatorManagerCrowd::EstimatorManagerCrowd(EstimatorManagerNew& em)
{
  // For now I'm going to try to refactor away the clone pattern only at the manager level.
  // i.e. not continue into the scalar_estimators and collectables
  main_estimator_ = UPtr<ScalarEstimatorBase>(em.main_estimator_->clone());
  for (const auto& est : em.scalar_ests_)
    scalar_estimators_.emplace_back(est->clone());
  for (const auto& upeb : em.operator_ests_)
    operator_ests_.emplace_back(upeb->spawnCrowdClone());
}

void EstimatorManagerCrowd::accumulate(const RefVector<MCPWalker>& walkers,
                                       const RefVector<ParticleSet>& psets,
                                       const RefVector<TrialWaveFunction>& wfns,
                                       RandomGenerator& rng)
{
  block_num_samples_ += walkers.size();
  for (MCPWalker& awalker : walkers)
    block_weight_ += awalker.Weight;
  main_estimator_->accumulate(walkers);
  int num_scalar_estimators = scalar_estimators_.size();
  for (int i = 0; i < num_scalar_estimators; ++i)
    scalar_estimators_[i]->accumulate(walkers);
  for (int i = 0; i < operator_ests_.size(); ++i)
    operator_ests_[i]->accumulate(walkers, psets, wfns, rng);
}


void EstimatorManagerCrowd::startBlock(int steps)
{
  for (auto& uope : operator_ests_)
    uope->startBlock(steps);
  block_num_samples_ = 0.0;
  block_weight_      = 0.0;
}

void EstimatorManagerCrowd::stopBlock()
{
  // the main estimator does more here.
}


} // namespace qmcplusplus
