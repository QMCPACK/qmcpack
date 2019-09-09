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
 *  TODO: Construct and initialize in thread execution space
 */
class Crowd
{
public:
  using MCPWalker = MCPopulation::MCPWalker;
  using GradType = QMCTraits::GradType;
  using RealType = QMCTraits::RealType;
  /** This is the data structure for walkers within a crowd
   */
  Crowd(EstimatorManagerBase& emb) : estimator_manager_crowd_(emb) {}

  /** With so many vectors allocate upfront.
   *
   *  could be premature optimization
   */
  void reserve(int crowd_size)
  {
    auto reserveCS = [crowd_size](auto& avector) {
                       avector.reserve(crowd_size); };
    reserveCS(mcp_walkers_);
    reserveCS(walker_elecs_);
    reserveCS(walker_twfs_);
    reserveCS(walker_hamiltonians_);
    reserveCS(grads_now_);
    reserveCS(grads_new_);
    reserveCS(ratios_);
  }
  
  void startRun() {}

  void startBlock(int steps) { estimator_manager_crowd_.startBlock(steps); }

  void addWalker(MCPWalker& walker, ParticleSet& elecs, TrialWaveFunction& twf, QMCHamiltonian& hamiltonian)
  {
    mcp_walkers_.push_back(walker);
    walker_elecs_.push_back(elecs);
    walker_twfs_.push_back(twf);
    walker_hamiltonians_.push_back(hamiltonian);
  };

  auto beginWalkers() { return mcp_walkers_.begin(); }
  auto endWalkers() { return mcp_walkers_.end(); }
  auto beginTrialWaveFunctions() { return walker_twfs_.begin(); }
  auto endTrialWaveFunctions() { return walker_twfs_.end(); }
  auto beginElectrons() { return walker_elecs_.begin(); }
  auto endElectrons() { return walker_elecs_.end(); }

  std::vector<std::reference_wrapper<ParticleSet>>& get_walker_elecs() { return walker_elecs_; }
  std::vector<std::reference_wrapper<TrialWaveFunction>>& get_walker_twfs() { return walker_twfs_; }
  std::vector<std::reference_wrapper<QMCHamiltonian>>& get_walker_hamiltonians() { return walker_hamiltonians_; }

  std::vector<GradType>& get_grads_now() { return grads_now_; }
  std::vector<GradType>& get_grads_new() { return grads_new_; }
  std::vector<TrialWaveFunction::PsiValueType>& get_ratios() { return ratios_; }
  
  int size() const { return mcp_walkers_.size(); }


  std::vector<RealType>& get_log_gf() { return log_gf_; }
  std::vector<RealType>& get_log_gb() { return log_gb_; }
  std::vector<RealType>& get_prob() { return prob_; }
  
  void start()
  {
    // These were cleared to 1.0 each loop by VMCUpdatePbyP advance walker
    // refactored code may depend on this initial value.
    std::fill(log_gf_.begin(), log_gf_.end(), 1.0);
    std::fill(log_gb_.begin(), log_gb_.end(), 1.0);
  }
  
private:
  std::vector<std::reference_wrapper<MCPWalker>> mcp_walkers_;
  std::vector<std::reference_wrapper<ParticleSet>> walker_elecs_;
  std::vector<std::reference_wrapper<TrialWaveFunction>> walker_twfs_;
  std::vector<std::reference_wrapper<QMCHamiltonian>> walker_hamiltonians_;
  /** @name Work Buffers
   *  @{
   *  There are many "local" variables in execution of the driver that are convenient to 
   *  place in a stl containers to use <alogorithm> style calls with.
   *  Eventually this will also us to allow some sort of appropriate parallel policy for these calls.
   */
  EstimatorManagerCrowd estimator_manager_crowd_;
  std::vector<GradType> grads_now_;
  std::vector<GradType> grads_new_;
  std::vector<TrialWaveFunction::PsiValueType> ratios_;
  std::vector<RealType> log_gf_;
  std::vector<RealType> log_gb_;
  std::vector<RealType> prob_;
  /** }@ */
};
} // namespace qmcplusplus
#endif
