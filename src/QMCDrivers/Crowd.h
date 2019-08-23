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
 */
class Crowd
{
public:
  using MCPWalker = MCPopulation::MCPWalker;
  /** This is the data structure for walkers within a crowd
   */
  struct Walkers
  {
  public:
    int live;
  };

  Crowd(EstimatorManagerBase emb) : estimator_manager_(emb) {}
  
  void startRun()
  {
    
  }

  void startBlock(int steps)
  {
    estimator_manager_.startBlock(steps);
  }

  void addWalker(MCPWalker& walker) { mcp_walkers_.push_back(walker); };
private:
  Walkers walkers_;
  std::vector<std::reference_wrapper<MCPWalker>> mcp_walkers_;
  
  EstimatorManagerCrowd estimator_manager_;
public:
};
} // namespace qmcplusplus
#endif
