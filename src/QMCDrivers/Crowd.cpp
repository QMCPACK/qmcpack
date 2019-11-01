//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/Crowd.h"

namespace qmcplusplus
{
void Crowd::clearWalkers()
{
  mcp_walkers_.clear();
  mcp_wfbuffers_.clear();
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
  mcp_wfbuffers_.push_back(walker.DataSet);
  walker_elecs_.push_back(elecs);
  walker_twfs_.push_back(twf);
  walker_hamiltonians_.push_back(hamiltonian);
};

void Crowd::loadWalkers()
{
  auto it_walker       = mcp_walkers_.begin();
  auto it_walker_elecs = walker_elecs_.begin();
  //flex walkers here
  while (it_walker != mcp_walkers_.end())
  {
    (*it_walker_elecs).get().loadWalker(*it_walker, true);
    ++it_walker;
    ++it_walker_elecs;
  }
}

void Crowd::startBlock(int num_steps)
{
  n_accept_ = 0;
  n_reject_ = 0;
  // VMCBatched does no nonlocal moves
  n_nonlocal_accept_ = 0;
  estimator_manager_crowd_.startBlock(num_steps);
}

} // namespace qmcplusplus
