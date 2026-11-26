//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include "Utilities/for_testing/Catch2Approx.h"

#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"

#include "type_traits/template_types.hpp"
#include "Particle/Walker.h"
#include "Estimators/EstimatorManagerNew.h"
#include "ParticleSetPool.h"
#include "WaveFunctionPool.h"
#include "HamiltonianPool.h"
#include "QMCDrivers/SFNBranch.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "QMCDrivers/tests/SetupPools.h"

namespace qmcplusplus
{
using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;
using RealType  = double;

namespace testing
{
class SetupSFNBranch
{
public:
  SetupSFNBranch(Communicate& comm) : comm_{comm} {}

  SetupSFNBranch() : comm_{*OHMMS::Controller} {}

  std::unique_ptr<SFNBranch> operator()(ParticleSet& pset, TrialWaveFunction& twf, QMCHamiltonian& ham)
  {
    pop_ = std::make_unique<MCPopulation>(1, comm_.rank(), pset, twf, ham);
    // MCPopulation owns it walkers it cannot just take refs so we just create and then update its walkers.
    pop_->createWalkers(2, walker_confs_);

    RefVector<MCPWalker> walkers = convertUPtrToRefVector(pop_->get_walkers());

    walkers[0].get().R[0] = 1.0;
    walkers[1].get().R[0] = 0.5;

    auto sfnb = std::make_unique<SFNBranch>(tau_, 1.0, DMCRefEnergyScheme::LIMITED_HISTORY);

    createMyNode(*sfnb, valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);

    sfnb->initParam(*pop_, 0, 0, false, false);

    return sfnb;
  }

private:
  void createMyNode(SFNBranch& sfnb, const char* xml)
  {
    doc_ = std::make_unique<Libxml2Document>();
    doc_->parseFromString(xml);
    sfnb.put(doc_->getRoot());
  }

  Communicate& comm_;
  UPtr<EstimatorManagerNew> emb_;
  WalkerConfigurations walker_confs_;
  UPtr<MCPopulation> pop_;
  UPtr<Libxml2Document> doc_;

  RealType tau_           = 1.0;
  int num_global_walkers_ = 16;
};

} // namespace testing

TEST_CASE("SFNBranch::branch(MCPopulation...)", "[drivers]")
{
  using namespace testing;
  SetupPools pools;
  SetupSFNBranch setup_sfnb(pools.comm);
  std::unique_ptr<SFNBranch> sfnb =
      setup_sfnb(*pools.particle_pool->getParticleSet("e"), pools.wavefunction_pool->getWaveFunction().value(),
                 pools.hamiltonian_pool->getHamiltonian().value());
}

TEST_CASE("SFNBranch fixed weight feedback and history", "[drivers][stochastic_comb]")
{
  testing::SetupPools pools;
  MCPopulation population(pools.comm.size(), pools.comm.rank(), *pools.particle_pool->getParticleSet("e"),
                          pools.wavefunction_pool->getWaveFunction().value(),
                          pools.hamiltonian_pool->getHamiltonian().value());
  WalkerConfigurations configs;
  population.createWalkers(2, configs);
  population.syncWalkersPerRank(pools.comm);
  const int count  = population.get_num_global_walkers();
  const int target = 3 * count / 2;
  for (auto scheme : {DMCRefEnergyScheme::LIMITED_HISTORY, DMCRefEnergyScheme::UNLIMITED_HISTORY})
    for (int warmup : {0, 2})
      for (bool fixed : {false, true})
      {
        SFNBranch branch(0.25, 0.4, scheme);
        Libxml2Document doc;
        const std::string xml = "<qmc><parameter name=\"warmupsteps\">" + std::to_string(warmup) +
            "</parameter><parameter name=\"energyUpdateInterval\">2</parameter>"
            "<parameter name=\"target_walkers\">" +
            std::to_string(target) + "</parameter><parameter name=\"warmupByReconfiguration\">" +
            (fixed ? "yes" : "no") + "</parameter></qmc>";
        REQUIRE(doc.parseFromString(xml));
        branch.put(doc.getRoot());
        branch.initParam(population, -1, 0.5, fixed, false);
        double expected_trial = -1;
        double sum_energy = 0, sum_weight = 0;
        for (int step = 0; step < warmup + 4; ++step)
        {
          MCDataType<double> ensemble;
          ensemble.Weight         = count * (step + 2.0);
          ensemble.Energy         = -2 - 0.1 * step;
          ensemble.Variance       = 0.5;
          ensemble.R2Accepted     = 1;
          ensemble.R2Proposed     = 1;
          ensemble.LivingFraction = 1;
          branch.updateParamAfterPopControl(ensemble, population.get_golden_electrons().getTotalNum());
          double expected_ref;
          if (step < warmup)
          {
            expected_ref   = ensemble.Energy;
            expected_trial = ensemble.Energy + std::log(double(target) / ensemble.Weight) / 0.25;
          }
          else
          {
            const double temporal_weight = scheme == DMCRefEnergyScheme::LIMITED_HISTORY ? ensemble.Weight : 1;
            sum_energy += temporal_weight * ensemble.Energy;
            sum_weight += temporal_weight;
            expected_ref = sum_energy / sum_weight;
            if ((step - warmup + 1) % 2 == 0)
              expected_trial = expected_ref + 0.4 * std::log(double(target) / ensemble.Weight);
          }
          CHECK(branch.getEref() == Approx(expected_ref));
          CHECK(branch.getEtrial() == Approx(expected_trial));
          CHECK(branch.getTauEff() == Approx(0.25));
          CHECK(population.get_num_global_walkers() == count);
        }
        MCDataType<double> invalid;
        invalid.Weight = 0;
        CHECK_THROWS_AS(branch.updateParamAfterPopControl(invalid, 2), UniformCommunicateError);
      }
}


} // namespace qmcplusplus
