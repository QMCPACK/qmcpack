//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch.hpp>

#include "Message/Communicate.h"

#include "type_traits/template_types.hpp"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/Walker.h"
#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/tests/FakeEstimator.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
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
class SetupSimpleFixedNodeBranch
{
public:
  SetupSimpleFixedNodeBranch(Communicate* comm)
  {
    comm_ = comm;
    emb_ = std::make_unique<EstimatorManagerBase>(comm_);
    emb_->add(std::make_unique<FakeEstimator>(), "fake");
  }

  SetupSimpleFixedNodeBranch()
  {
    comm_ = OHMMS::Controller;
    emb_ = std::make_unique<EstimatorManagerBase>(comm_);
    emb_->add(std::make_unique<FakeEstimator>(), "fake");
  }

  SimpleFixedNodeBranch operator()()
  {
    const SimulationCell simulation_cell;
    mcwc_ = std::make_unique<MCWalkerConfiguration>(simulation_cell);
    mcwc_->setName("electrons");

    mcwc_->create(1);
    mcwc_->R[0][0] = 0.0;
    mcwc_->R[0][1] = 1.0;
    mcwc_->R[0][2] = 2.0;

    MCPWalker w1(1);
    w1.R[0] = 1.0;
    MCPWalker w2(1);
    w2.R[0] = 0.5;

    mcwc_->createWalkers(2);
    (*mcwc_)[0]->R = w1.R;
    (*mcwc_)[1]->R = w2.R;

    SimpleFixedNodeBranch sfnb(tau_, num_global_walkers_);

    sfnb.setEstimatorManager(std::move(emb_));

    createMyNode(sfnb, valid_dmc_input_sections[valid_dmc_input_dmc_index]);

    // WalkerController is created as a side effect.
    sfnb.initWalkerController(*mcwc_, false, false);

    sfnb.checkParameters(*mcwc_);

    sfnb.branch(0, *mcwc_);

    return sfnb;
  }

  
private:

  void createMyNode(SimpleFixedNodeBranch& sfnb, const char* xml)
  {
    doc_ = std::make_unique<Libxml2Document>();
    doc_->parseFromString(xml);
    sfnb.myNode = doc_->getRoot();
  }
  
  Communicate* comm_;
  UPtr<EstimatorManagerBase> emb_;
  UPtr<MCPopulation> pop_;
  UPtr<Libxml2Document> doc_;
  
  RealType tau_           = 1.0;
  int num_global_walkers_ = 16;

  //Legacy
  UPtr<MCWalkerConfiguration> mcwc_;
};

}

TEST_CASE("SimpleFixedNodeBranch::branch(MCWC...)", "[drivers][legacy]")
{
  using namespace testing;
  SetupSimpleFixedNodeBranch setup_sfnb;
  SimpleFixedNodeBranch sfnb = setup_sfnb();  
  
  // \todo: check walker ID's here. might need more walkers to make this significant.

}

} // namespace qmcplusplus
