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
#include "Particle/MCPopulation.h"
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
  /** This is the data structure for walkers within a crowd
   */
  Crowd(EstimatorManagerBase& emb) : estimator_manager_crowd_(emb) {}

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
  
  int size() const { return mcp_walkers_.size(); }

private:
  std::vector<std::reference_wrapper<MCPWalker>> mcp_walkers_;
  std::vector<std::reference_wrapper<ParticleSet>> walker_elecs_;
  std::vector<std::reference_wrapper<TrialWaveFunction>> walker_twfs_;
  std::vector<std::reference_wrapper<QMCHamiltonian>> walker_hamiltonians_;
  EstimatorManagerCrowd estimator_manager_crowd_;

public:
};
} // namespace qmcplusplus
#endif
