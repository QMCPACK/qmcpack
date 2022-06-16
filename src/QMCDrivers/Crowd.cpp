//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "Crowd.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{

Crowd::Crowd(EstimatorManagerNew& emb, const MultiWalkerDispatchers& dispatchers)
    : dispatchers_(dispatchers), estimator_manager_crowd_(emb)
{}

Crowd::Crowd(EstimatorManagerNew& emb,
             ParticleSet& pset_electrons,
             TrialWaveFunction& twf,
             QMCHamiltonian& ham,
             const MultiWalkerDispatchers& dispatchers)
    : dispatchers_(dispatchers),
      driverwalker_resource_collection_(createCrowdResources(pset_electrons, twf, ham)),
      estimator_manager_crowd_(emb)
{}

Crowd::~Crowd() = default;

void Crowd::clearWalkers()
{
  // We're clearing the refs to the objects not the referred to objects.
  mcp_walkers_.clear();
  walker_elecs_.clear();
  walker_twfs_.clear();
  walker_hamiltonians_.clear();
}

DriverWalkerResourceCollection Crowd::createCrowdResources(ParticleSet& pset,
                                                           TrialWaveFunction& twf,
                                                           QMCHamiltonian& ham)
{
  DriverWalkerResourceCollection resources;
  pset.createResource(resources.pset_res);
  twf.createResource(resources.twf_res);
  ham.createResource(resources.ham_res);
  return resources;
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

void Crowd::setRNGForHamiltonian(RandomGenerator& rng)
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
