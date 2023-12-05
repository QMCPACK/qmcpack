//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "Crowd.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
Crowd::Crowd(EstimatorManagerNew& emb,
             const DriverWalkerResourceCollection& driverwalker_res,
             const ParticleSet& pset,
             const TrialWaveFunction& twf,
             const QMCHamiltonian& ham,
             const MultiWalkerDispatchers& dispatchers)
  : dispatchers_(dispatchers), driverwalker_resource_collection_(driverwalker_res), estimator_manager_crowd_(emb)
{
  if (emb.areThereListeners())
  {
    // By creating a tempory QMCHamiltonian before walkers are distributed to the crowds we insure that
    // after construction that each crowd has a valid driverwalker_resource_collection_.ham_res.
    // We can't use the golden hamiltonian to do this because QMCHamiltonian and mw_res_ are 1 to 1 or 1 to 0.
    //
    // Violating the incapsulation of the QMChamiltonian mw_res_ could make this very efficent but doesn't seem
    // necessary since this happens num crowds times per section.
    //
    // QMCHamiltonian makes quite a smell with its non const pset and twf constructor
    // arguments.
    ParticleSet pset_temp(pset);
    UPtr<TrialWaveFunction> twf_temp(twf.makeClone(pset_temp));
    UPtr<QMCHamiltonian> ham_temp(ham.makeClone(pset_temp, *twf_temp));
    RefVectorWithLeader<QMCHamiltonian> ham_list{*ham_temp, {*ham_temp}};
    ResourceCollectionTeamLock<QMCHamiltonian> ham_lock(driverwalker_resource_collection_.ham_res, ham_list);
    estimator_manager_crowd_.registerListeners(ham_list);
  }
}

Crowd::~Crowd() = default;

void Crowd::clearWalkers()
{
  // We're clearing the refs to the objects not the referred to objects.
  mcp_walkers_.clear();
  walker_elecs_.clear();
  walker_twfs_.clear();
  walker_hamiltonians_.clear();
}

void Crowd::reserve(int crowd_size)
{
  auto reserveCS = [crowd_size](auto& avector) { avector.reserve(crowd_size); };
  reserveCS(mcp_walkers_);
  reserveCS(walker_elecs_);
  reserveCS(walker_twfs_);
  reserveCS(walker_hamiltonians_);
}

void Crowd::addWalker(MCPWalker& walker, ParticleSet& elecs, TrialWaveFunction& twf, QMCHamiltonian& hamiltonian)
{
  mcp_walkers_.push_back(walker);
  walker_elecs_.push_back(elecs);
  walker_twfs_.push_back(twf);
  walker_hamiltonians_.push_back(hamiltonian);
};

void Crowd::setRNGForHamiltonian(RandomBase<FullPrecRealType>& rng)
{
  for (QMCHamiltonian& ham : walker_hamiltonians_)
    ham.setRandomGenerator(&rng);
}

void Crowd::startBlock(int num_steps)
{
  n_accept_ = 0;
  n_reject_ = 0;
  // VMCBatched does no nonlocal moves
  n_nonlocal_accept_ = 0;
  estimator_manager_crowd_.startBlock(num_steps);
}

void Crowd::stopBlock() { estimator_manager_crowd_.stopBlock(); }

} // namespace qmcplusplus
