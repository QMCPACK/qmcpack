//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
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
class EstimatorManagerNew;
class EstimatorManagerCrowd;

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
  EnergyDensityTest(Communicate* comm,
                    int num_walkers,
                    std::function<MockGoldWalkerElements(Communicate*, RuntimeOptions)> make_gold_elem,
                    bool generate_test_data = false);

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

  std::vector<ParticleSet::ParticlePos> deterministic_rs_2 = {{{-0.2708307179, -0.7751078697, -1.507662187},
                                                               {-1.767935983, -2.267908723, 0.802928773},
                                                               {2.638861876, 3.189295837, -0.1583547099},
                                                               {1.967142627, 3.168934617, 0.1348337978},
                                                               {-0.3154878117, -0.8555576282, -0.615753412},
                                                               {-0.5064132498, -0.04847732172, -1.201739871},
                                                               {-0.01394719094, 1.845506534, 1.376975874},
                                                               {-0.3800458355, 2.255976232, 1.629132391}},
                                                              {{-2.489025257, 0.8482124128, 1.434512321},
                                                               {-0.0767998772, -0.2892376071, -0.152022511},
                                                               {-1.041033132, 1.93429606, 2.190503856},
                                                               {4.682283313, 3.679345099, 3.037770313},
                                                               {1.122318726, -1.044722727, 0.2493639639},
                                                               {-0.1709426256, 0.91126951, 0.0234699242},
                                                               {3.128878401, 0.7625744583, 1.585080796},
                                                               {2.594682586, 1.887794866, 2.210192703}},
                                                              {{0.5225379708, 1.072364434, -2.89798479},
                                                               {-0.3527898799, 0.376540703, 0.3485486555},
                                                               {2.59226874, 1.335303561, 1.280898488},
                                                               {-0.3819189956, 0.5358460986, 2.720153363},
                                                               {-0.2665835197, 0.5974395087, -0.5135086651},
                                                               {-0.6419828646, 0.4538051387, -1.057858504},
                                                               {1.420183152, 1.676011905, -0.01909935602},
                                                               {2.175727247, 1.327399178, 0.9107993414}},
                                                              {{-0.3681009669, 1.271148301, 0.7084059291},
                                                               {-1.596724176, -0.03391919965, -0.7795750272},
                                                               {2.935133797, 2.369577977, 2.555840556},
                                                               {0.3823055891, 2.169392989, 0.3343588146},
                                                               {-0.989526036, 0.4286373805, -1.75878267},
                                                               {0.614791086, -1.391095985, -1.367703815},
                                                               {2.957505805, 0.3219368438, 2.851501195},
                                                               {1.826007474, 2.738471999, 0.5792334772}}};

  std::vector<ParticleSet::ParticlePos> deterministic_rs_3 = {{{1.1740665, 1.708365287, -0.7221926758},
                                                               {0.7191311468, 0.8161897766, 0.361988661},
                                                               {2.055987137, 0.8035338848, 1.735882921},
                                                               {1.62884783, 0.9373587674, 2.687607336},
                                                               {0.3155039079, -0.6311718493, 1.072053495},
                                                               {-1.562583557, -1.362662151, -0.2322816531},
                                                               {0.3864303135, 2.609428725, 1.75511768},
                                                               {0.4422101671, 1.213026597, 0.5986030429}},
                                                              {{0.8923285479, 0.9345148512, 0.4521097489},
                                                               {0.8909405863, 0.2229965821, -1.282158307},
                                                               {1.258704658, 2.484271561, 1.010671276},
                                                               {2.522248942, 3.671887703, 1.783684111},
                                                               {-0.7675185706, -1.893068888, -2.247128921},
                                                               {-2.644468122, -0.209719895, 0.8985644648},
                                                               {1.862680294, 3.282609522, 1.490447074},
                                                               {2.798469359, 1.10809932, 3.481222147}},
                                                              {{1.600027786, -0.9474347311, 0.4711368114},
                                                               {-0.7611051151, 0.5765779072, 0.1967857514},
                                                               {2.136350635, 3.188981458, 1.4786544},
                                                               {1.406956883, 2.237787843, 1.404264641},
                                                               {0.7537326229, 0.01526880603, 1.846934965},
                                                               {0.7467097515, -0.758435228, 0.3651868226},
                                                               {2.312927579, 0.7089260252, 0.6434838012},
                                                               {2.505633299, 1.490758679, 2.607601639}},
                                                              {{0.7726521842, 0.3962054602, 0.3567443323},
                                                               {-1.338373662, 1.704015782, -0.7761974804},
                                                               {2.167978672, 2.34190624, 1.220030335},
                                                               {1.77832023, 1.308655557, 1.265443042},
                                                               {-2.017466086, -1.691870446, 0.4042196747},
                                                               {0.1987181218, 0.4657786182, 1.286565222},
                                                               {1.718174661, 3.822324813, 0.8313790978},
                                                               {1.338128839, 1.399575664, 1.921516492}}};

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
