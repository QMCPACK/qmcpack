//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "Configuration.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/DeepQMC/DeepQMCBridge.h"
#include "QMCWaveFunctions/DeepQMC/DeepQMCWF.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/TWFdispatcher.h"
#include "QMCWaveFunctions/TWFGrads.hpp"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "Utilities/RuntimeOptions.h"
#include "ResourceCollection.h"

namespace qmcplusplus
{
using RealType = DeepQMCBridge::RealType;

namespace
{
class RecordingDeepQMCBridge : public DeepQMCBridge
{
public:
  BatchResult evaluateLogBatch(const std::vector<RealType>& ion_coords,
                               const std::vector<RealType>& electron_coords,
                               int mol_idx,
                               int batch_size,
                               int n_elec) const override
  {
    call_count++;
    last_ion_coords      = ion_coords;
    last_electron_coords = electron_coords;
    last_mol_idx         = mol_idx;
    last_batch_size      = batch_size;
    last_n_elec          = n_elec;

    BatchResult result;
    result.log_values.resize(batch_size);
    result.grad_log_values.resize(static_cast<std::size_t>(batch_size) * n_elec * OHMMS_DIM);
    result.lap_log_values.resize(static_cast<std::size_t>(batch_size) * n_elec);

    for (int iw = 0; iw < batch_size; ++iw)
    {
      result.log_values[iw] = logValue(electron_coords, iw, n_elec);
      for (int iat = 0; iat < n_elec; ++iat)
      {
        const std::size_t particle_offset = static_cast<std::size_t>(iw) * n_elec + iat;
        for (int d = 0; d < OHMMS_DIM; ++d)
          result.grad_log_values[particle_offset * OHMMS_DIM + d] = 100.0 * iw + 10.0 * iat + d;
        result.lap_log_values[particle_offset] = 1000.0 * iw + iat;
      }
    }
    return result;
  }

  virtual RealType logValue(const std::vector<RealType>& electron_coords, int iw, int n_elec) const
  {
    return 10.0 + iw;
  }

  mutable int call_count = 0;
  mutable std::vector<RealType> last_ion_coords;
  mutable std::vector<RealType> last_electron_coords;
  mutable int last_mol_idx    = -1;
  mutable int last_batch_size = -1;
  mutable int last_n_elec     = -1;
};

class CoordinateLogDeepQMCBridge : public RecordingDeepQMCBridge
{
public:
  RealType logValue(const std::vector<RealType>& electron_coords, int iw, int n_elec) const override
  {
    RealType log_value      = 0.0;
    const std::size_t begin = static_cast<std::size_t>(iw) * n_elec * OHMMS_DIM;
    for (int i = 0; i < n_elec * OHMMS_DIM; ++i)
      log_value += electron_coords[begin + i];
    return log_value;
  }
};

ParticleSet makeElectrons(const SimulationCell& simulation_cell,
                          std::initializer_list<ParticleSet::SingleParticlePos> positions)
{
  ParticleSet electrons(simulation_cell);
  electrons.setName("e");
  electrons.create({static_cast<int>(positions.size())});
  int iat = 0;
  for (const auto& pos : positions)
    electrons.R[iat++] = pos;
  electrons.update();
  return electrons;
}

ParticleSet makeIons(const SimulationCell& simulation_cell)
{
  ParticleSet ions(simulation_cell);
  ions.setName("ion0");
  ions.create({1});
  ions.R[0] = {1.0, 2.0, 3.0};
  ions.update();
  return ions;
}

class UnsupportedAllParticleWF : public WaveFunctionComponent
{
public:
  std::string getClassName() const override { return "UnsupportedAllParticleWF"; }
  LogValue evaluateLog(const ParticleSet&, ParticleSet::ParticleGradient&, ParticleSet::ParticleLaplacian&) override
  {
    return LogValue(0);
  }
  void acceptMove(ParticleSet&, int, bool) override {}
  void restore(int) override {}
  PsiValue ratio(ParticleSet&, int) override { return PsiValue(1); }
  void registerData(ParticleSet&, WFBufferType&) override {}
  LogValue updateBuffer(ParticleSet&, WFBufferType&, bool) override { return LogValue(0); }
  void copyFromBuffer(ParticleSet&, WFBufferType&) override {}
  void evaluateDerivatives(ParticleSet&, const OptVariables&, Vector<ValueType>&, Vector<ValueType>&) override {}
};
} // namespace

TEST_CASE("DeepQMCWF batched evaluateLog", "[wavefunction][deepqmc]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions  = makeIons(simulation_cell);
  ParticleSet elec0 = makeElectrons(simulation_cell, {{0.0, 0.1, 0.2}, {1.0, 1.1, 1.2}});
  ParticleSet elec1 = makeElectrons(simulation_cell, {{2.0, 2.1, 2.2}, {3.0, 3.1, 3.2}});

  auto bridge_owner = std::make_unique<RecordingDeepQMCBridge>();
  auto* bridge      = bridge_owner.get();
  DeepQMCWF comp0("deep", ions, std::move(bridge_owner), 7);
  DeepQMCWF comp1("deep", ions, std::make_unique<RecordingDeepQMCBridge>(), 7);

  ParticleSet::ParticleGradient G0, G1;
  ParticleSet::ParticleLaplacian L0, L1;
  G0.resize(elec0.getTotalNum());
  G1.resize(elec1.getTotalNum());
  L0.resize(elec0.getTotalNum());
  L1.resize(elec1.getTotalNum());
  G0 = 0.0;
  G1 = 0.0;
  L0 = 0.0;
  L1 = 0.0;

  RefVectorWithLeader<WaveFunctionComponent> wfc_list(comp0);
  wfc_list.push_back(comp0);
  wfc_list.push_back(comp1);
  RefVectorWithLeader<ParticleSet> p_list(elec0);
  p_list.push_back(elec0);
  p_list.push_back(elec1);
  RefVector<ParticleSet::ParticleGradient> G_list;
  G_list.push_back(G0);
  G_list.push_back(G1);
  RefVector<ParticleSet::ParticleLaplacian> L_list;
  L_list.push_back(L0);
  L_list.push_back(L1);

  comp0.mw_evaluateLog(wfc_list, p_list, G_list, L_list);

  CHECK(bridge->call_count == 1);
  CHECK(bridge->last_mol_idx == 7);
  CHECK(bridge->last_batch_size == 2);
  CHECK(bridge->last_n_elec == 2);
  REQUIRE(bridge->last_ion_coords.size() == 3);
  CHECK(bridge->last_ion_coords[0] == Approx(1.0));
  CHECK(bridge->last_ion_coords[1] == Approx(2.0));
  CHECK(bridge->last_ion_coords[2] == Approx(3.0));

  const std::vector<RealType> expected_electron_coords{0.0, 0.1, 0.2, 1.0, 1.1, 1.2, 2.0, 2.1, 2.2, 3.0, 3.1, 3.2};
  REQUIRE(bridge->last_electron_coords.size() == expected_electron_coords.size());
  for (int i = 0; i < expected_electron_coords.size(); ++i)
    CHECK(bridge->last_electron_coords[i] == Approx(expected_electron_coords[i]));

  CHECK(std::real(comp0.get_log_value()) == Approx(10.0));
  CHECK(std::real(comp1.get_log_value()) == Approx(11.0));

  CHECK(G0[0][0] == Approx(0.0));
  CHECK(G0[0][1] == Approx(1.0));
  CHECK(G0[0][2] == Approx(2.0));
  CHECK(G0[1][0] == Approx(10.0));
  CHECK(G0[1][1] == Approx(11.0));
  CHECK(G0[1][2] == Approx(12.0));
  CHECK(L0[0] == Approx(0.0));
  CHECK(L0[1] == Approx(1.0));

  CHECK(G1[0][0] == Approx(100.0));
  CHECK(G1[0][1] == Approx(101.0));
  CHECK(G1[0][2] == Approx(102.0));
  CHECK(G1[1][0] == Approx(110.0));
  CHECK(G1[1][1] == Approx(111.0));
  CHECK(G1[1][2] == Approx(112.0));
  CHECK(L1[0] == Approx(1000.0));
  CHECK(L1[1] == Approx(1001.0));
}

TEST_CASE("DeepQMCWF all-particle masked evaluation scatters into the full crowd", "[wavefunction][deepqmc]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions  = makeIons(simulation_cell);
  ParticleSet elec0 = makeElectrons(simulation_cell, {{0.0, 0.0, 0.0}});
  ParticleSet elec1 = makeElectrons(simulation_cell, {{1.0, 0.0, 0.0}});
  ParticleSet elec2 = makeElectrons(simulation_cell, {{2.0, 0.0, 0.0}});
  auto bridge_owner = std::make_unique<RecordingDeepQMCBridge>();
  auto* bridge      = bridge_owner.get();
  DeepQMCWF comp0("deep", ions, std::move(bridge_owner), 7);
  DeepQMCWF comp1("deep", ions, std::make_unique<RecordingDeepQMCBridge>(), 7);
  DeepQMCWF comp2("deep", ions, std::make_unique<RecordingDeepQMCBridge>(), 7);
  ParticleSet::ParticleGradient G0(1), G1(1), G2(1);
  ParticleSet::ParticleLaplacian L0(1), L1(1), L2(1);
  G0 = 0.0;
  G1 = 17.0;
  G2 = 0.0;
  L0 = 0.0;
  L1 = 19.0;
  L2 = 0.0;
  RefVectorWithLeader<WaveFunctionComponent> wfc_list(comp0, {comp0, comp1, comp2});
  RefVectorWithLeader<ParticleSet> p_list(elec0, {elec0, elec1, elec2});
  RefVector<ParticleSet::ParticleGradient> G_list;
  G_list.push_back(G0);
  G_list.push_back(G1);
  G_list.push_back(G2);
  RefVector<ParticleSet::ParticleLaplacian> L_list;
  L_list.push_back(L0);
  L_list.push_back(L1);
  L_list.push_back(L2);

  comp0.mw_evaluateLogAllParticles(wfc_list, p_list, G_list, L_list, {true, false, true});

  CHECK(bridge->call_count == 1);
  CHECK(bridge->last_batch_size == 2);
  CHECK(bridge->last_electron_coords[0] == Approx(0.0));
  CHECK(bridge->last_electron_coords[3] == Approx(2.0));
  CHECK(std::real(comp0.get_log_value()) == Approx(10.0));
  CHECK(std::real(comp1.get_log_value()) == Approx(0.0));
  CHECK(std::real(comp2.get_log_value()) == Approx(11.0));
  CHECK(G0[0][0] == Approx(0.0));
  CHECK(G1[0][0] == Approx(17.0));
  CHECK(G2[0][0] == Approx(100.0));
  CHECK(L1[0] == Approx(19.0));
  CHECK(L2[0] == Approx(1000.0));
}

TEST_CASE("all-particle TWF dispatch recovers rejected walkers only", "[wavefunction][deepqmc]")
{
  const SimulationCell simulation_cell;
  RuntimeOptions runtime_options;
  ParticleSet ions  = makeIons(simulation_cell);
  ParticleSet elec0 = makeElectrons(simulation_cell, {{0.0, 0.0, 0.0}});
  ParticleSet elec1 = makeElectrons(simulation_cell, {{2.0, 0.0, 0.0}});
  auto bridge_owner = std::make_unique<CoordinateLogDeepQMCBridge>();
  auto* bridge      = bridge_owner.get();
  TrialWaveFunction twf0(runtime_options, "deep0");
  twf0.addComponent(std::make_unique<DeepQMCWF>("DNN", ions, std::move(bridge_owner), 7));
  TrialWaveFunction twf1(runtime_options, "deep1");
  twf1.addComponent(std::make_unique<DeepQMCWF>("DNN", ions, std::make_unique<CoordinateLogDeepQMCBridge>(), 7));
  RefVectorWithLeader<TrialWaveFunction> wf_list(twf0, {twf0, twf1});
  RefVectorWithLeader<ParticleSet> p_list(elec0, {elec0, elec1});
  ResourceCollection resources("deepqmc_all_particle_resources");
  elec0.createResource(resources);
  ResourceCollectionTeamLock<ParticleSet> lock(resources, p_list);
  MCCoords<CoordsType::POS> displacements(2);
  displacements.positions = {{1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
  std::vector<bool> valid(2);
  ParticleSet::mw_makeMoveAllParticles(p_list, displacements, valid);
  TWFdispatcher dispatcher(true);
  dispatcher.flex_evaluateLogAllParticles(wf_list, p_list, {true, true});
  const auto accepted_log = twf0.getLogPsi();
  const ParticleSet::ParticleGradient accepted_g(elec0.G);
  const ParticleSet::ParticleLaplacian accepted_l(elec0.L);
  ParticleSet::mw_accept_rejectMoveAllParticles(p_list, {true, false});
  dispatcher.flex_accept_rejectMoveAllParticles(wf_list, p_list, {true, false});
  CHECK(bridge->call_count == 2);
  CHECK(bridge->last_batch_size == 1);
  CHECK(twf0.getLogPsi() == Approx(accepted_log));
  CHECK(twf1.getLogPsi() == Approx(2.0));
  for (int idim = 0; idim < OHMMS_DIM; ++idim)
    CHECK(elec0.G[0][idim] == Approx(accepted_g[0][idim]));
  CHECK(elec0.L[0] == Approx(accepted_l[0]));

  ParticleSet::mw_makeMoveAllParticles(p_list, displacements, valid);
  dispatcher.flex_evaluateLogAllParticles(wf_list, p_list, {true, true});
  const int calls_before_accept = bridge->call_count;
  ParticleSet::mw_accept_rejectMoveAllParticles(p_list, {true, true});
  dispatcher.flex_accept_rejectMoveAllParticles(wf_list, p_list, {true, true});
  CHECK(bridge->call_count == calls_before_accept);
}

TEST_CASE("all-particle WFC capability and TWF dispatcher are gated", "[wavefunction][deepqmc]")
{
  const SimulationCell simulation_cell;
  RuntimeOptions runtime_options;
  ParticleSet elec = makeElectrons(simulation_cell, {{0.0, 0.0, 0.0}});
  TrialWaveFunction twf(runtime_options, "unsupported");
  twf.addComponent(std::make_unique<UnsupportedAllParticleWF>());
  RefVectorWithLeader<TrialWaveFunction> wf_list(twf, {twf});
  RefVectorWithLeader<ParticleSet> p_list(elec, {elec});
  REQUIRE_THROWS_AS(TrialWaveFunction::mw_evaluateLogAllParticles(wf_list, p_list, {true}), std::runtime_error);
  REQUIRE_THROWS_AS(TrialWaveFunction::mw_accept_rejectMoveAllParticles(wf_list, p_list, {false}), std::runtime_error);
  TWFdispatcher dispatcher(false);
  REQUIRE_THROWS_AS(dispatcher.flex_evaluateLogAllParticles(wf_list, p_list, {true}), std::runtime_error);
  REQUIRE_THROWS_AS(dispatcher.flex_accept_rejectMoveAllParticles(wf_list, p_list, {true}), std::runtime_error);
}

TEST_CASE("DeepQMCWF batched PbyP methods", "[wavefunction][deepqmc]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions  = makeIons(simulation_cell);
  ParticleSet elec0 = makeElectrons(simulation_cell, {{0.0, 0.1, 0.2}, {1.0, 1.1, 1.2}});
  ParticleSet elec1 = makeElectrons(simulation_cell, {{2.0, 2.1, 2.2}, {3.0, 3.1, 3.2}});

  auto bridge_owner = std::make_unique<CoordinateLogDeepQMCBridge>();
  auto* bridge      = bridge_owner.get();
  DeepQMCWF comp0("deep", ions, std::move(bridge_owner), 7);
  DeepQMCWF comp1("deep", ions, std::make_unique<CoordinateLogDeepQMCBridge>(), 7);

  ParticleSet::ParticleGradient G0(elec0.getTotalNum()), G1(elec1.getTotalNum());
  ParticleSet::ParticleLaplacian L0(elec0.getTotalNum()), L1(elec1.getTotalNum());
  G0 = 0.0;
  G1 = 0.0;
  L0 = 0.0;
  L1 = 0.0;
  comp0.evaluateLog(elec0, G0, L0);
  comp1.evaluateLog(elec1, G1, L1);
  const auto old_log0 = std::real(comp0.get_log_value());
  const auto old_log1 = std::real(comp1.get_log_value());

  RefVectorWithLeader<WaveFunctionComponent> wfc_list(comp0);
  wfc_list.push_back(comp0);
  wfc_list.push_back(comp1);
  RefVectorWithLeader<ParticleSet> p_list(elec0);
  p_list.push_back(elec0);
  p_list.push_back(elec1);

  std::vector<WaveFunctionComponent::GradType> grads(2);
  comp0.mw_evalGrad(wfc_list, p_list, 1, grads);
  CHECK(bridge->last_batch_size == 2);
  CHECK(bridge->last_electron_coords[0] == Approx(0.0));
  CHECK(bridge->last_electron_coords[6] == Approx(2.0));
  CHECK(grads[0][0] == Approx(10.0));
  CHECK(grads[0][2] == Approx(12.0));
  CHECK(grads[1][0] == Approx(110.0));
  CHECK(grads[1][2] == Approx(112.0));

  elec0.makeMove(1, ParticleSet::SingleParticlePos{0.5, 0.0, 0.0});
  elec1.makeMove(1, ParticleSet::SingleParticlePos{0.0, 0.25, 0.0});

  std::vector<WaveFunctionComponent::PsiValue> ratios;
  std::vector<WaveFunctionComponent::GradType> grad_new(2);
  comp0.mw_ratioGrad(wfc_list, p_list, 1, ratios, grad_new);

  const RealType new_log0 = old_log0 + 0.5;
  const RealType new_log1 = old_log1 + 0.25;
  CHECK(bridge->last_batch_size == 2);
  const std::vector<RealType> expected_active_coords{0.0, 0.1, 0.2, 1.5, 1.1, 1.2, 2.0, 2.1, 2.2, 3.0, 3.35, 3.2};
  REQUIRE(bridge->last_electron_coords.size() == expected_active_coords.size());
  for (int i = 0; i < expected_active_coords.size(); ++i)
    CHECK(bridge->last_electron_coords[i] == Approx(expected_active_coords[i]));
  CHECK(ratios[0] == Approx(std::exp(new_log0 - old_log0)));
  CHECK(ratios[1] == Approx(std::exp(new_log1 - old_log1)));
  CHECK(grad_new[0][1] == Approx(11.0));
  CHECK(grad_new[1][1] == Approx(111.0));

  comp0.mw_accept_rejectMove(wfc_list, p_list, 1, {true, false}, true);
  CHECK(std::real(comp0.get_log_value()) == Approx(new_log0));
  CHECK(std::real(comp1.get_log_value()) == Approx(old_log1));

  comp0.mw_calcRatio(wfc_list, p_list, 1, ratios);
  CHECK(ratios[0] == Approx(1.0));
  CHECK(ratios[1] == Approx(std::exp(new_log1 - old_log1)));
}

TEST_CASE("TrialWaveFunction dispatches DeepQMC batched PbyP methods", "[wavefunction][deepqmc]")
{
  const SimulationCell simulation_cell;
  RuntimeOptions runtime_options;
  ParticleSet ions  = makeIons(simulation_cell);
  ParticleSet elec0 = makeElectrons(simulation_cell, {{0.0, 0.1, 0.2}, {1.0, 1.1, 1.2}});
  ParticleSet elec1 = makeElectrons(simulation_cell, {{2.0, 2.1, 2.2}, {3.0, 3.1, 3.2}});

  auto bridge_owner = std::make_unique<CoordinateLogDeepQMCBridge>();
  auto* bridge      = bridge_owner.get();
  TrialWaveFunction twf0(runtime_options, "deepqmc0");
  twf0.addComponent(std::make_unique<DeepQMCWF>("DNN", ions, std::move(bridge_owner), 7));
  TrialWaveFunction twf1(runtime_options, "deepqmc1");
  twf1.addComponent(std::make_unique<DeepQMCWF>("DNN", ions, std::make_unique<CoordinateLogDeepQMCBridge>(), 7));

  RefVectorWithLeader<TrialWaveFunction> wf_list(twf0);
  wf_list.push_back(twf0);
  wf_list.push_back(twf1);
  RefVectorWithLeader<ParticleSet> p_list(elec0);
  p_list.push_back(elec0);
  p_list.push_back(elec1);

  TrialWaveFunction::mw_evaluateLog(wf_list, p_list);
  const auto old_log0 = twf0.getLogPsi();
  const auto old_log1 = twf1.getLogPsi();

  TWFGrads<CoordsType::POS> grads_now(2);
  TrialWaveFunction::mw_evalGrad(wf_list, p_list, 1, grads_now);
  CHECK(bridge->last_batch_size == 2);
  CHECK(grads_now.grads_positions[0][2] == Approx(12.0));
  CHECK(grads_now.grads_positions[1][2] == Approx(112.0));

  elec0.makeMove(1, ParticleSet::SingleParticlePos{0.5, 0.0, 0.0});
  elec1.makeMove(1, ParticleSet::SingleParticlePos{0.0, 0.25, 0.0});

  std::vector<TrialWaveFunction::PsiValue> ratios;
  TWFGrads<CoordsType::POS> grads_new(2);
  TrialWaveFunction::mw_calcRatioGrad(wf_list, p_list, 1, ratios, grads_new);

  const RealType new_log0 = old_log0 + 0.5;
  const RealType new_log1 = old_log1 + 0.25;
  CHECK(bridge->last_batch_size == 2);
  CHECK(ratios[0] == Approx(std::exp(new_log0 - old_log0)));
  CHECK(ratios[1] == Approx(std::exp(new_log1 - old_log1)));
  CHECK(grads_new.grads_positions[0][1] == Approx(11.0));
  CHECK(grads_new.grads_positions[1][1] == Approx(111.0));

  TrialWaveFunction::mw_accept_rejectMove(wf_list, p_list, 1, {true, false}, true);
  CHECK(twf0.getLogPsi() == Approx(new_log0));
  CHECK(twf1.getLogPsi() == Approx(old_log1));
}

TEST_CASE("DeepQMCWF batched evaluateGL delegates to one batch", "[wavefunction][deepqmc]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions  = makeIons(simulation_cell);
  ParticleSet elec0 = makeElectrons(simulation_cell, {{0.0, 0.1, 0.2}, {1.0, 1.1, 1.2}});
  ParticleSet elec1 = makeElectrons(simulation_cell, {{2.0, 2.1, 2.2}, {3.0, 3.1, 3.2}});

  auto bridge_owner = std::make_unique<RecordingDeepQMCBridge>();
  auto* bridge      = bridge_owner.get();
  DeepQMCWF comp0("deep", ions, std::move(bridge_owner), 7);
  DeepQMCWF comp1("deep", ions, std::make_unique<RecordingDeepQMCBridge>(), 7);

  ParticleSet::ParticleGradient G0(elec0.getTotalNum()), G1(elec1.getTotalNum());
  ParticleSet::ParticleLaplacian L0(elec0.getTotalNum()), L1(elec1.getTotalNum());
  G0 = 0.0;
  G1 = 0.0;
  L0 = 0.0;
  L1 = 0.0;

  RefVectorWithLeader<WaveFunctionComponent> wfc_list(comp0);
  wfc_list.push_back(comp0);
  wfc_list.push_back(comp1);
  RefVectorWithLeader<ParticleSet> p_list(elec0);
  p_list.push_back(elec0);
  p_list.push_back(elec1);
  RefVector<ParticleSet::ParticleGradient> G_list;
  G_list.push_back(G0);
  G_list.push_back(G1);
  RefVector<ParticleSet::ParticleLaplacian> L_list;
  L_list.push_back(L0);
  L_list.push_back(L1);

  comp0.mw_evaluateGL(wfc_list, p_list, G_list, L_list, true);

  CHECK(bridge->call_count == 1);
  CHECK(bridge->last_batch_size == 2);
  CHECK(std::real(comp0.get_log_value()) == Approx(10.0));
  CHECK(std::real(comp1.get_log_value()) == Approx(11.0));
  CHECK(G1[1][2] == Approx(112.0));
  CHECK(L1[1] == Approx(1001.0));
}

TEST_CASE("PythonDeepQMCBridge calls Python batch bridge", "[wavefunction][deepqmc]")
{
  namespace fs              = std::filesystem;
  const fs::path bridge_dir = fs::temp_directory_path() / "qmcpack_deepqmc_bridge_test";
  fs::remove_all(bridge_dir);
  fs::create_directories(bridge_dir);

  std::ofstream bridge_py(bridge_dir / "deepqmc_infer_bridge.py");
  bridge_py << R"PY(
class DeepQMCInferBridge:
    def __init__(self, model_path):
        self.model_path = model_path

    def compute_log_gl(self, nuclear_coords, electron_coords, mol_idx, batch_size, n_elec):
        assert self.model_path == 'fake-model'
        assert nuclear_coords == [1.0, 2.0, 3.0]
        assert electron_coords == [0.0, 0.1, 0.2, 1.0, 1.1, 1.2]
        assert batch_size == 1
        assert n_elec == 2
        return [20.0 + mol_idx], [0.0, 1.0, 2.0, 10.0, 11.0, 12.0], [100.0, 101.0]
)PY";
  bridge_py.close();

  auto bridge       = makePythonDeepQMCBridge("fake-model", bridge_dir.string());
  const auto result = bridge->evaluateLogBatch({1.0, 2.0, 3.0}, {0.0, 0.1, 0.2, 1.0, 1.1, 1.2}, 5, 1, 2);

  REQUIRE(result.log_values.size() == 1);
  REQUIRE(result.grad_log_values.size() == 6);
  REQUIRE(result.lap_log_values.size() == 2);
  CHECK(result.log_values[0] == Approx(25.0));
  CHECK(result.grad_log_values[5] == Approx(12.0));
  CHECK(result.lap_log_values[1] == Approx(101.0));

  fs::remove_all(bridge_dir);
}

TEST_CASE("WaveFunctionFactory builds DeepQMC component from XML", "[wavefunction][deepqmc]")
{
  namespace fs              = std::filesystem;
  const fs::path bridge_dir = fs::temp_directory_path() / "qmcpack_deepqmc_factory_test";
  fs::remove_all(bridge_dir);
  fs::create_directories(bridge_dir);

  std::ofstream bridge_py(bridge_dir / "deepqmc_infer_bridge.py");
  bridge_py << R"PY(
class DeepQMCInferBridge:
    def __init__(self, model_path):
        self.model_path = model_path

    def compute_log_gl(self, nuclear_coords, electron_coords, mol_idx, batch_size, n_elec):
        assert self.model_path == 'factory-model'
        assert nuclear_coords == [0.0, 0.0, 0.0]
        assert electron_coords == [0.0, 0.1, 0.2, 1.0, 1.1, 1.2,
                                   2.0, 2.1, 2.2, 3.0, 3.1, 3.2]
        assert mol_idx == 9
        assert batch_size == 2
        assert n_elec == 2
        return [30.0, 31.0], [0.0, 1.0, 2.0, 10.0, 11.0, 12.0,
                              100.0, 101.0, 102.0, 110.0, 111.0, 112.0], [200.0, 201.0, 300.0, 301.0]
)PY";
  bridge_py.close();

  const SimulationCell simulation_cell;
  WaveFunctionFactory::PSetMap particle_set_map;

  auto ions = std::make_unique<ParticleSet>(simulation_cell);
  ions->setName("ion0");
  ions->create({1});
  ions->R[0] = {0.0, 0.0, 0.0};
  ions->update();
  particle_set_map.emplace("ion0", std::move(ions));

  auto elec0 = std::make_unique<ParticleSet>(simulation_cell);
  elec0->setName("e");
  elec0->create({2});
  elec0->R[0] = {0.0, 0.1, 0.2};
  elec0->R[1] = {1.0, 1.1, 1.2};
  elec0->update();
  particle_set_map.emplace("e", std::move(elec0));

  ParticleSet elec1 = makeElectrons(simulation_cell, {{2.0, 2.1, 2.2}, {3.0, 3.1, 3.2}});

  std::string wavefunction_xml = "<wavefunction name=\"psi0\" target=\"e\">"
                                 "<deepqmc name=\"DNN\" source=\"ion0\" model=\"factory-model\" mol_idx=\"9\" "
                                 "python_module_path=\"" +
      bridge_dir.string() +
      "\"/>"
      "</wavefunction>";
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(wavefunction_xml.c_str()));

  WaveFunctionFactory wff(*particle_set_map["e"], particle_set_map, OHMMS::Controller);
  RuntimeOptions runtime_options;
  auto twf0 = wff.buildTWF(doc.getRoot(), runtime_options);
  REQUIRE(twf0 != nullptr);
  REQUIRE(twf0->size() == 1);
  auto twf1 = twf0->makeClone(elec1);

  RefVectorWithLeader<TrialWaveFunction> wf_list(*twf0);
  wf_list.push_back(*twf0);
  wf_list.push_back(*twf1);
  RefVectorWithLeader<ParticleSet> p_list(*particle_set_map["e"]);
  p_list.push_back(*particle_set_map["e"]);
  p_list.push_back(elec1);

  TrialWaveFunction::mw_evaluateLog(wf_list, p_list);

  CHECK(twf0->getLogPsi() == Approx(30.0));
  CHECK(twf1->getLogPsi() == Approx(31.0));
  CHECK(twf0->G[1][2] == Approx(12.0));
  CHECK(twf1->G[0][0] == Approx(100.0));
  CHECK(twf0->L[1] == Approx(201.0));
  CHECK(twf1->L[1] == Approx(301.0));

  fs::remove_all(bridge_dir);
}

TEST_CASE("PythonDeepQMCBridge can load a real DeepQMC He checkpoint", "[wavefunction][deepqmc]")
{
  const char* checkpoint  = std::getenv("DEEPQMC_HE_CHECKPOINT");
  const char* python_path = std::getenv("DEEPQMC_PYTHON_SITE_PACKAGES");
  if (checkpoint == nullptr || python_path == nullptr)
  {
    SUCCEED("Set DEEPQMC_HE_CHECKPOINT and DEEPQMC_PYTHON_SITE_PACKAGES to run the optional real DeepQMC bridge test");
    return;
  }

  auto bridge       = makePythonDeepQMCBridge(checkpoint, python_path);
  const auto result = bridge->evaluateLogBatch({0.0, 0.0, 0.0}, {0.1, 0.0, 0.0, -0.1, 0.0, 0.0}, 0, 1, 2);

  REQUIRE(result.log_values.size() == 1);
  REQUIRE(result.grad_log_values.size() == 6);
  REQUIRE(result.lap_log_values.size() == 2);
  CHECK(std::abs(result.log_values[0]) < 1.0e100);
  for (const auto value : result.grad_log_values)
    CHECK(std::abs(value) < 1.0e100);
  for (const auto value : result.lap_log_values)
    CHECK(std::abs(value) < 1.0e100);
}

TEST_CASE("DeepQMCWF single walker delegates to batched evaluateLog", "[wavefunction][deepqmc]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions = makeIons(simulation_cell);
  ParticleSet elec = makeElectrons(simulation_cell, {{0.0, 0.1, 0.2}, {1.0, 1.1, 1.2}});

  auto bridge_owner = std::make_unique<RecordingDeepQMCBridge>();
  auto* bridge      = bridge_owner.get();
  DeepQMCWF comp("deep", ions, std::move(bridge_owner), 3);

  ParticleSet::ParticleGradient G;
  ParticleSet::ParticleLaplacian L;
  G.resize(elec.getTotalNum());
  L.resize(elec.getTotalNum());
  G = 0.0;
  L = 0.0;

  auto log_value = comp.evaluateLog(elec, G, L);

  CHECK(bridge->call_count == 1);
  CHECK(bridge->last_batch_size == 1);
  CHECK(bridge->last_n_elec == 2);
  CHECK(bridge->last_mol_idx == 3);
  CHECK(std::real(log_value) == Approx(10.0));
  CHECK(G[1][2] == Approx(12.0));
  CHECK(L[1] == Approx(1.0));
}

} // namespace qmcplusplus
