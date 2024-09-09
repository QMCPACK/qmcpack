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

  /// Canned test data so rng differences don't cause fails.
  std::vector<ParticleSet::ParticlePos> deterministic_rs_ = {
      {
          {0.515677886, 0.9486072745, -1.17423246},
          {-0.3166678423, 1.376550506, 1.445290031},
          {1.96071365, 2.47265689, 1.051449486},
          {0.745853269, 0.5551359072, 4.080774681},
          {-0.3515016103, -0.5192222523, 0.9941510909},
          {-0.8354426872, 0.7071638258, -0.3409843552},
          {0.4386044751, 1.237378731, 2.331874152},
          {2.125850717, 0.3221067321, 0.5825731561},
      },
      {
          {-0.4633736785, 0.06318772224, -0.8580153742},
          {-1.174926354, -0.6276503679, 0.07458759314},
          {1.327618206, 2.085829379, 1.415749862},
          {0.9114727103, 0.1789183931, -0.08135540251},
          {-2.267908723, 0.802928773, 0.9522812957},
          {1.502715257, -1.84493529, 0.2805620469},
          {3.168934617, 0.1348337978, 1.371092768},
          {0.8310229518, 1.070827168, 1.18016733},
      },
      {
          {-0.04847732172, -1.201739871, -1.700527771},
          {0.1589259538, -0.3096047065, -2.066626415},
          {2.255976232, 1.629132391, -0.8024446773},
          {2.534792993, 3.121092901, 1.609780703},
          {-0.2892376071, -0.152022511, -2.727613712},
          {0.2477154804, 0.5039232765, 2.995702733},
          {3.679345099, 3.037770313, 2.808899306},
          {0.6418578532, 1.935944544, 1.515637954},
      },
      {
          {0.91126951, 0.0234699242, 1.442297821},
          {-0.9240061217, -0.1014997844, 0.9081020061},
          {1.887794866, 2.210192703, 2.209118551},
          {2.758945014, -1.21140421, 1.3337907},
          {0.376540703, 0.3485486555, 0.9056881595},
          {-0.3512770187, -0.4056820917, -2.068499576},
          {0.5358460986, 2.720153363, 1.41999706},
          {2.284020089, 1.173071915, 1.044597715},
      },
  };
};

} // namespace testing
} // namespace qmcplusplus

#endif
