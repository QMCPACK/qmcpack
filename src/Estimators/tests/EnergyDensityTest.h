//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ENERGYDENSITYTEST_H
#define QMCPLUSPLUS_ENERGYDENSITYTEST_H

#include "EnergyDensityEstimator.h"
#include "ResourceCollection.h"
#include "Utilities/ProjectData.h"
#include "MockGoldWalkerElements.h"

namespace qmcplusplus
{
namespace testing
{

/** This class is here to allow sharing of the code to stand up an evironment to test
 *  the energy density estimator.
 */
class EnergyDensityTest
{
public:
  using MCPWalker = typename OperatorEstBase::MCPWalker;

  EnergyDensityTest(Communicate* comm, int num_walkers, bool generate_test_data = false);
  
  RefVector<ParticleSet> getPsetRefs() { return makeRefVector<ParticleSet>(psets_); }
  NEEnergyDensityEstimator& getEnergyDensityEstimator();

  MockGoldWalkerElements& getGoldElements() { return gold_elem_; }
  
  RefVectorWithLeader<MCPWalker> getWalkerList();
  RefVectorWithLeader<ParticleSet> getPSetList();
  RefVectorWithLeader<QMCHamiltonian> getHamList();
  RefVectorWithLeader<TrialWaveFunction> getTwfList();

  ResourceCollection& getPSetRes();
  ResourceCollection& getHamRes();
  ResourceCollection& getTwfRes();
  
  template <typename LIST>
  ResourceCollectionTeamLock<LIST> makeTeamLock(ResourceCollection& res, RefVectorWithLeader<LIST> list) {
    return ResourceCollectionTeamLock<LIST>(res, list);
  }
private:
  ProjectData test_project_;
  MockGoldWalkerElements gold_elem_;

  UPtrVector<QMCHamiltonian> hams_;
  UPtrVector<TrialWaveFunction> twfs_;
  std::vector<ParticleSet> psets_;
  std::vector<MCPWalker> walkers_;

  ResourceCollection pset_res_{"test_pset_res"};
  ResourceCollection ham_res_{"test_ham_res"};
  ResourceCollection wfc_res_{"test_wfc_res"};

  UPtr<NEEnergyDensityEstimator> eden_est_;
  UPtr<EnergyDensityInput> edein_;
};

} // namespace testing
} // namespace qmcplusplus

#endif
