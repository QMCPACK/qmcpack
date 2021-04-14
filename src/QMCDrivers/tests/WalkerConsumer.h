//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_WALKERCONSUMER_H
#define QMCPLUSPLUS_WALKERCONSUMER_H

#include <vector>

#include "Particle/Walker.h"

namespace qmcplusplus
{
class ResourceCollection;

namespace testing
{
/** mock class to avoid testing dependency between Crowd and MCPopulation
 *
 *  Also example of minimum client of MCPopulation::redistributeWalkers
 */
class WalkerConsumer
{
public:
  std::vector<std::reference_wrapper<Walker<QMCTraits, PtclOnLatticeTraits>>> walkers;
  std::vector<std::reference_wrapper<ParticleSet>> walker_elecs_;
  std::vector<std::reference_wrapper<TrialWaveFunction>> walker_twfs_;
  std::vector<std::reference_wrapper<QMCHamiltonian>> walker_hamiltonians_;

  void initializeResources(const ResourceCollection& twf_resource) {}

  void addWalker(Walker<QMCTraits, PtclOnLatticeTraits>& walker,
                 ParticleSet& elecs,
                 TrialWaveFunction& twf,
                 QMCHamiltonian& hamiltonian)
  {
    walkers.push_back(walker);
    walker_elecs_.push_back(elecs);
    walker_twfs_.push_back(twf);
    walker_hamiltonians_.push_back(hamiltonian);
  }

  void clearWalkers()
  {
    // We're clearing the refs to the objects not the referred to objects.
    walkers.clear();
    walker_elecs_.clear();
    walker_twfs_.clear();
    walker_hamiltonians_.clear();
  }
};

} // namespace testing
} // namespace qmcplusplus
#endif
