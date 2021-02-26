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
#include "ResourceCollection.h"

namespace qmcplusplus
{
Crowd::Crowd(EstimatorManagerNew& emb, const MultiWalkerDispatchers& dispatchers)
    : dispatchers_(dispatchers), estimator_manager_crowd_(emb)
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

void Crowd::loadWalkers()
{
  for (int i = 0; i < mcp_walkers_.size(); ++i)
    walker_elecs_[i].get().loadWalker(mcp_walkers_[i], true);
}

void Crowd::setRNGForHamiltonian(RandomGenerator_t& rng)
{
  for (QMCHamiltonian& ham : walker_hamiltonians_)
    ham.setRandomGenerator(&rng);
}

void Crowd::initializeResources(const ResourceCollection& twf_resource)
{
  if (!twfs_shared_resource_)
    twfs_shared_resource_ = std::make_unique<ResourceCollection>(twf_resource);
}

void Crowd::lendResources(size_t receiver)
{
  if (walker_twfs_.size() > 0)
  {
    if (twfs_shared_resource_->is_lent())
      throw std::runtime_error("Crowd::lendResources resources are out already!");
    if (receiver >= walker_twfs_.size())
      throw std::runtime_error("Crowd::takebackResources receiver out of bound!");
    twfs_shared_resource_->lock();
    twfs_shared_resource_->rewind();
    walker_twfs_[receiver].get().acquireResource(*twfs_shared_resource_);
  }
}

void Crowd::takebackResources(size_t receiver)
{
  if (walker_twfs_.size() > 0)
  {
    if (!twfs_shared_resource_->is_lent())
      throw std::runtime_error("Crowd::takebackResources resources was not lent out!");
    if (receiver >= walker_twfs_.size())
      throw std::runtime_error("Crowd::takebackResources receiver out of bound!");
    twfs_shared_resource_->rewind();
    walker_twfs_[receiver].get().releaseResource(*twfs_shared_resource_);
    twfs_shared_resource_->unlock();
  }
}

void Crowd::startBlock(int num_steps)
{
  if (this->size() == 0)
    return;
  n_accept_ = 0;
  n_reject_ = 0;
  // VMCBatched does no nonlocal moves
  n_nonlocal_accept_ = 0;
  estimator_manager_crowd_.startBlock(num_steps);
}

void Crowd::stopBlock()
{
  if (this->size() == 0)
    return;
  estimator_manager_crowd_.stopBlock();
}

} // namespace qmcplusplus
