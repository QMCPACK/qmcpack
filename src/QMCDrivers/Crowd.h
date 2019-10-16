//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from VMC.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CROWD_H
#define QMCPLUSPLUS_CROWD_H

#include <vector>
#include "QMCDrivers/MCPopulation.h"
#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/EstimatorManagerCrowd.h"

namespace qmcplusplus
{
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
  using WFBuffer         = MCPopulation::WFBuffer;
  using GradType         = QMCTraits::GradType;
  using RealType         = QMCTraits::RealType;
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  /** This is the data structure for walkers within a crowd
   */
  Crowd(EstimatorManagerBase& emb) : estimator_manager_crowd_(emb) {}

  /** Because so many vectors allocate them upfront.
   *
   *  could be premature optimization
   */
  void reserve(int crowd_size);
  void startRun() {}
  void startBlock(int steps);

  EstimatorManagerCrowd& get_estimator_manager_crowd() { return estimator_manager_crowd_; }
  void addWalker(MCPWalker& walker, ParticleSet& elecs, TrialWaveFunction& twf, QMCHamiltonian& hamiltonian);

  void loadWalkers();

  /** clears log_gf and log_gb
   *
   *  This is a legacy "call"
   *  TODO: likely remove log_gf and log_gb from crowd
   *        seems unlikely these should be persistent state.
   */
  void clearResults();
  
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
    estimator_manager_crowd_.accumulate(global_walkers, mcp_walkers_, walker_elecs_);
  }

  auto beginWalkers() { return mcp_walkers_.begin(); }
  auto endWalkers() { return mcp_walkers_.end(); }
  auto beginTrialWaveFunctions() { return walker_twfs_.begin(); }
  auto endTrialWaveFunctions() { return walker_twfs_.end(); }
  auto beginElectrons() { return walker_elecs_.begin(); }
  auto endElectrons() { return walker_elecs_.end(); }

  RefVector<MCPWalker>& get_walkers() { return mcp_walkers_; }
  std::vector<std::reference_wrapper<ParticleSet>>& get_walker_elecs() { return walker_elecs_; }
  std::vector<std::reference_wrapper<TrialWaveFunction>>& get_walker_twfs() { return walker_twfs_; }
  std::vector<std::reference_wrapper<QMCHamiltonian>>& get_walker_hamiltonians() { return walker_hamiltonians_; }

  std::vector<GradType>& get_grads_now() { return grads_now_; }
  std::vector<GradType>& get_grads_new() { return grads_new_; }
  std::vector<TrialWaveFunction::PsiValueType>& get_ratios() { return ratios_; }
  std::vector<RealType>& get_log_gf() { return log_gf_; }
  std::vector<RealType>& get_log_gb() { return log_gb_; }
  std::vector<RealType>& get_prob() { return prob_; }
  RefVector<WFBuffer>& get_mcp_wfbuffers() { return mcp_wfbuffers_; }
  const EstimatorManagerCrowd& get_estimator_manager_crowd() const { return estimator_manager_crowd_; }
  int size() const { return mcp_walkers_.size(); }

  void incReject() { ++n_reject_; }
  void incAccept() { ++n_accept_; }
  void incNonlocalAccept(int n = 1) { n_nonlocal_accept_ += n; }
  FullPrecRealType get_accept_ratio() const
  {
    return [](FullPrecRealType accept, FullPrecRealType reject) -> FullPrecRealType {
      return accept / (accept + reject);
    }(n_accept_, n_reject_);
  }
  unsigned long get_nonlocal_accept() { return n_nonlocal_accept_; }
private:
  void resizeResults(int crowd_size);
  /** @name Walker Vectors
   *
   *  A single index into these ordered lists constitue a complete 
   *  walker context.
   * @{
   */
  RefVector<MCPWalker> mcp_walkers_;
  RefVector<WFBuffer> mcp_wfbuffers_;
  RefVector<ParticleSet> walker_elecs_;
  RefVector<TrialWaveFunction> walker_twfs_;
  RefVector<QMCHamiltonian> walker_hamiltonians_;
  /** }@ */
  
  EstimatorManagerCrowd estimator_manager_crowd_;
  /** @name Work Buffers
   *  @{
   *  There are many "local" variables in execution of the driver that are convenient to 
   *  place in a stl containers to use <alogorithm> style calls with.
   *  Eventually this will also us to allow some sort of appropriate parallel policy for these calls.
   */
  std::vector<GradType> grads_now_;
  std::vector<GradType> grads_new_;
  std::vector<TrialWaveFunction::PsiValueType> ratios_;
  std::vector<RealType> log_gf_;
  std::vector<RealType> log_gb_;
  std::vector<RealType> prob_;
  /** }@ */

  /** @name Step State
   * 
   *  Should be per walker? 
   *  @{
   */
  unsigned long n_reject_ = 0;
  unsigned long n_accept_ = 0;
  unsigned long n_nonlocal_accept_ = 0;
  /** @} */
};
} // namespace qmcplusplus
#endif
