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

#include <catch.hpp>

#include "Message/Communicate.h"

#include "type_traits/template_types.hpp"
#include "Particle/Walker.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Estimators/tests/FakeEstimator.h"
#include "QMCDrivers/SFNBranch.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/Crowd.h"
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
  SetupSFNBranch(Communicate* comm)
  {
    comm_                   = comm;
    emb_                    = std::make_unique<EstimatorManagerNew>(comm_);
    FakeEstimator* fake_est = new FakeEstimator;
    emb_->add(fake_est, "fake");
  }

  SetupSFNBranch()
  {
    comm_                   = OHMMS::Controller;
    emb_                    = std::make_unique<EstimatorManagerNew>(comm_);
    FakeEstimator* fake_est = new FakeEstimator;
    emb_->add(fake_est, "fake");
  }

  SFNBranch operator()(ParticleSet& p_set, TrialWaveFunction& twf, QMCHamiltonian& ham)
  {
    pop_ = std::make_unique<MCPopulation>(1, comm_->rank(), walker_confs, &p_set, &twf, &ham);
    // MCPopulation owns it walkers it cannot just take refs so we just create and then update its walkers.
    pop_->createWalkers(2);

    RefVector<MCPWalker> walkers = convertUPtrToRefVector(pop_->get_walkers());

    walkers[0].get().R[0] = 1.0;
    walkers[1].get().R[0] = 0.5;

    SFNBranch sfnb(tau_, num_global_walkers_);

    sfnb.setEstimatorManager(emb_.get());

    createMyNode(sfnb, valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);

    sfnb.initWalkerController(*pop_, false, false);


    sfnb.checkParameters(pop_->get_num_global_walkers(), walkers);

    UPtrVector<Crowd> crowds;
    crowds.emplace_back(std::make_unique<Crowd>(*emb_));
    crowds.emplace_back(std::make_unique<Crowd>(*emb_));

    sfnb.branch(0, *pop_);

    return sfnb;
  }

private:
  void createMyNode(SFNBranch& sfnb, const char* xml)
  {
    doc_ = std::make_unique<Libxml2Document>();
    doc_->parseFromString(xml);
    sfnb.myNode = doc_->getRoot();
  }

  Communicate* comm_;
  UPtr<EstimatorManagerNew> emb_;
  WalkerConfigurations walker_confs;
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
  SFNBranch sfnb = setup_sfnb(*pools.particle_pool->getParticleSet("e"), *pools.wavefunction_pool->getPrimary(),
                              *pools.hamiltonian_pool->getPrimary());
}


} // namespace qmcplusplus
