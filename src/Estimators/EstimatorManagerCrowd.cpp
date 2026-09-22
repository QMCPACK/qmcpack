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
  main_estimator_ = UPtr<ScalarEstimatorBase>(em.main_estimator_->clone());
  for (const auto& est : em.scalar_ests_)
    scalar_estimators_.emplace_back(est->clone());
  for (const auto& upeb : em.operator_ests_)
    operator_ests_.emplace_back(upeb->spawnCrowdClone());
}

void EstimatorManagerCrowd::accumulate(const RefVector<MCPWalker>& walkers,
                                       const RefVector<ParticleSet>& psets,
                                       const RefVector<TrialWaveFunction>& wfns,
                                       const RefVector<QMCHamiltonian>& hams,
                                       RandomBase<FullPrecRealType>& rng)
{
  block_num_samples_ += walkers.size();
  for (MCPWalker& awalker : walkers)
    block_weight_ += awalker.Weight;
  main_estimator_->accumulate(walkers);
  int num_scalar_estimators = scalar_estimators_.size();
  for (int i = 0; i < num_scalar_estimators; ++i)
    scalar_estimators_[i]->accumulate(walkers);
  for (int i = 0; i < operator_ests_.size(); ++i)
    operator_ests_[i]->accumulate(walkers, psets, wfns, hams, rng);
}

void EstimatorManagerCrowd::registerListeners(const RefVectorWithLeader<QMCHamiltonian>& ham_list)
{
  for (auto& estimator : operator_ests_)
    if (estimator->isListenerRequired())
      estimator->registerListeners(ham_list.getLeader());
}

void EstimatorManagerCrowd::startBlock(int steps)
{
  for (auto& uope : operator_ests_)
    uope->startBlock(steps);
  block_num_samples_ = 0.0;
  block_weight_      = 0.0;
  vmc_previous_weight_ = 0.0;
  vmc_data_.clear();
  vmc_data_.reserve(steps);
}

void EstimatorManagerCrowd::stopBlock()
{
  for (auto& uope : operator_ests_)
    uope->stopBlock();
}

void EstimatorManagerCrowd::recordVMCStep(unsigned long accepted, unsigned long rejected)
{
  std::vector<RealType> row;
  auto append_and_clear = [&row](ScalarEstimatorBase& estimator) {
    for (auto& scalar : estimator.scalars)
    {
      row.push_back(scalar.result());
      scalar.clear();
    }
  };
  append_and_clear(*main_estimator_);
  for (auto& estimator : scalar_estimators_)
    append_and_clear(*estimator);
  row.push_back(block_weight_ - vmc_previous_weight_);
  row.push_back(static_cast<RealType>(accepted));
  row.push_back(static_cast<RealType>(rejected));
  vmc_previous_weight_ = block_weight_;
  vmc_data_.emplace_back(std::move(row));
}


} // namespace qmcplusplus
