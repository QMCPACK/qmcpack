//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CROWD_H
#define QMCPLUSPLUS_CROWD_H

#include <vector>
#include "QMCDrivers/MCPopulation.h"
#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/EstimatorManagerCrowd.h"
#include "RandomGenerator.h"
#include "MultiWalkerDispatchers.h"

namespace qmcplusplus
{
// forward declaration
class ResourceCollection;

/** Driver synchronized step context
 * 
 *  assumed to live inside the drivers scope
 *  Represents the walker contexts exlusively operated on by
 *  one concurrent "worker" context from a total population
 *  owned by a MCPopulation 
 *
 *  TODO: Maybe construct and initialize in thread execution space
 *        @ye doubts this is important if rank is confined to one NUMA
 */
class Crowd
{
public:
  using MCPWalker        = MCPopulation::MCPWalker;
  using GradType         = QMCTraits::GradType;
  using RealType         = QMCTraits::RealType;
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  /** This is the data structure for walkers within a crowd
   */
  Crowd(EstimatorManagerNew& emb, const MultiWalkerDispatchers& dispatchers);

  ~Crowd();
  /** Because so many vectors allocate them upfront.
   *
   *  could be premature optimization
   */
  void reserve(int crowd_size);
  void startRun() {}
  void startBlock(int steps);
  void stopBlock();

  EstimatorManagerCrowd& get_estimator_manager_crowd() { return estimator_manager_crowd_; }
  void addWalker(MCPWalker& walker, ParticleSet& elecs, TrialWaveFunction& twf, QMCHamiltonian& hamiltonian);

  void loadWalkers();

  /** Clears all walker vectors
   *
   *  Unless you are _redistributing_ walkers to crowds don't
   *  call this. Prefer allowing the crowd lifecycle to roll.
   *  DO NOT move constructor code into here and call from constructor
   *  That is legacy "reset" pattern is considered deprecated
   */
  void clearWalkers();

  void accumulate(int global_walkers)
  {
    if (this->size() == 0)
      return;
    estimator_manager_crowd_.accumulate(global_walkers, mcp_walkers_, walker_elecs_);
  }

  void setRNGForHamiltonian(RandomGenerator_t& rng);

  /** initialize crowd-owned resources shared by walkers in the crowd
   */
  void initializeResources(const ResourceCollection& twf_resource);
  /** lend crowd-owned resources to the crowd leader walker
   * Note: use RAII CrowdResourceLock whenever possible
   */
  void lendResources(size_t receiver);
  /** take back crowd-owned resources from the crowd leader walker
   * Note: use RAII CrowdResourceLock whenever possible
   */
  void takebackResources(size_t receiver);

  auto beginWalkers() { return mcp_walkers_.begin(); }
  auto endWalkers() { return mcp_walkers_.end(); }
  auto beginTrialWaveFunctions() { return walker_twfs_.begin(); }
  auto endTrialWaveFunctions() { return walker_twfs_.end(); }
  auto beginElectrons() { return walker_elecs_.begin(); }
  auto endElectrons() { return walker_elecs_.end(); }

  const RefVector<MCPWalker>& get_walkers() const { return mcp_walkers_; }
  const RefVector<ParticleSet>& get_walker_elecs() const { return walker_elecs_; }
  const RefVector<TrialWaveFunction>& get_walker_twfs() const { return walker_twfs_; }
  const RefVector<QMCHamiltonian>& get_walker_hamiltonians() const { return walker_hamiltonians_; }

  const EstimatorManagerCrowd& get_estimator_manager_crowd() const { return estimator_manager_crowd_; }

  ResourceCollection& getTWFSharedResource() { return *twfs_shared_resource_; }

  int size() const { return mcp_walkers_.size(); }

  void incReject() { ++n_reject_; }
  void incAccept() { ++n_accept_; }
  void incNonlocalAccept(int n = 1) { n_nonlocal_accept_ += n; }
  unsigned long get_nonlocal_accept() { return n_nonlocal_accept_; }
  unsigned long get_accept() { return n_accept_; }
  unsigned long get_reject() { return n_reject_; }

  const MultiWalkerDispatchers& dispatchers_;

private:
  /** @name Walker Vectors
   *
   *  A single index into these ordered lists constitutes a complete 
   *  walker context.
   * @{
   */
  RefVector<MCPWalker> mcp_walkers_;
  RefVector<ParticleSet> walker_elecs_;
  RefVector<TrialWaveFunction> walker_twfs_;
  RefVector<QMCHamiltonian> walker_hamiltonians_;
  /** }@ */

  EstimatorManagerCrowd estimator_manager_crowd_;

  std::unique_ptr<ResourceCollection> twfs_shared_resource_;

  /** @name Step State
   * 
   *  Should be per walker? 
   *  @{
   */
  unsigned long n_reject_          = 0;
  unsigned long n_accept_          = 0;
  unsigned long n_nonlocal_accept_ = 0;
  /** @} */
};

/** Lock for a crowd lending and taking back shared resource to its consumer objects.
 */
class CrowdResourceLock
{
public:
  CrowdResourceLock(Crowd& locked_crowd, size_t receiver = 0) : locked_crowd_(locked_crowd), receiver_(receiver)
  {
    locked_crowd_.lendResources(receiver_);
  }

  ~CrowdResourceLock() { locked_crowd_.takebackResources(receiver_); }

  CrowdResourceLock(const CrowdResourceLock&) = delete;
  CrowdResourceLock(CrowdResourceLock&&)      = delete;

private:
  Crowd& locked_crowd_;
  const size_t receiver_;
};

} // namespace qmcplusplus
#endif
