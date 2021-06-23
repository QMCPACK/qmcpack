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
  SetupSFNBranch(Communicate* comm)
  {
    comm_ = comm;
    emb_  = std::make_unique<EstimatorManagerNew>(comm_);
  }

  SetupSFNBranch()
  {
    comm_ = OHMMS::Controller;
    emb_  = std::make_unique<EstimatorManagerNew>(comm_);
  }

  std::unique_ptr<SFNBranch> operator()(ParticleSet& p_set, TrialWaveFunction& twf, QMCHamiltonian& ham)
  {
    pop_ = std::make_unique<MCPopulation>(1, comm_->rank(), walker_confs_, &p_set, &twf, &ham);
    // MCPopulation owns it walkers it cannot just take refs so we just create and then update its walkers.
    pop_->createWalkers(2);

    RefVector<MCPWalker> walkers = convertUPtrToRefVector(pop_->get_walkers());

    walkers[0].get().R[0] = 1.0;
    walkers[1].get().R[0] = 0.5;

    auto sfnb = std::make_unique<SFNBranch>(tau_, num_global_walkers_);

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

  Communicate* comm_;
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
      setup_sfnb(*pools.particle_pool->getParticleSet("e"), *pools.wavefunction_pool->getPrimary(),
                 *pools.hamiltonian_pool->getPrimary());
}


} // namespace qmcplusplus
