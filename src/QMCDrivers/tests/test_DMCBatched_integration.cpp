//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "RuntimeOptions.h"
#include "catch.hpp"
#include <filesystem>

#include "QMCDrivers/QMCDriverFactory.h"
#include "DMC/DMCBatched.h"
#include "ValidQMCInputSections.h"
#include <MockGoldWalkerElements.h>
#include <EstimatorInputDelegates.h>
#include "DMCBatchedTestAccessor.h"
#include <RandomNumberControl.h>
#include <MakeEstimatorManager.h>
#include "DMC/DMCContextForSteps.h"
#include "reference_prefix.h"

namespace qmcplusplus
{

constexpr bool generate_new_seed  = false;
constexpr bool generate_test_data = true;

void copyRefFile(std::string_view ref_file, std::string_view ref_dest)
{
  std::filesystem::copy_file(ref_file, ref_dest, std::filesystem::copy_options::overwrite_existing);
}

TEST_CASE("DMCBatched::estimator_measurement_period", "[drivers]")
{
  using namespace testing;
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;
  outputManager.pause();
  using DMCInput = qmcplusplus::testing::DmcInputs;
  Libxml2Document doc;
  bool okay = doc.parseFromString(DMCInput::getXml(DMCInput::valid::EMPERIOD));
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();

  RuntimeOptions run_time_options;

  auto mgwe = makeGoldWalkerElementsWithEEEI(comm, run_time_options);

  QMCDriverFactory driver_factory(test_project);
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
  std::unique_ptr<QMCDriverInterface> qmc_driver;
  auto* qmc_system = mgwe.particle_pool.getWalkerSet("e");

  // This code repeats code in driver factory because it really only
  // exists because of legacy and erases driver type.
  QMCDriverInput qmcdriver_input;
  DMCDriverInput dmcdriver_input;
  qmcdriver_input.readXML(node);
  dmcdriver_input.readXML(node);

  if constexpr (generate_new_seed)
  {
    // This may not work for parallel hdf5
    RandomNumberControl::write(RandomNumberControl::getChildrenRefs(), "seeds_for_DMCBatched", comm);
    if (comm->rank() == 0)
    {
      std::string rand_file{"seeds_for_DMCBatched.random.h5"};
      std::string random_dest_file{"rank_count_" + std::to_string(comm->size()) + "_" + rand_file};
      copyRefFile(rand_file, random_dest_file);
    }
  }
  else
  {
    std::string rand_file{"seeds_for_DMCBatched"};
    std::string random_dest_file{"rank_count_" + std::to_string(comm->size()) + "_" + rand_file};
    RandomNumberControl::getChildrenRefs();
    RandomNumberControl::read(random_dest_file, comm);
  }

  auto random_child_refs = RandomNumberControl::getChildrenRefs();
  DMCBatched dmc_batched(test_project, std::move(qmcdriver_input),
                         makeEstimatorManager(std::nullopt, qmcdriver_input.get_estimator_manager_input(),
                                              mgwe.pset_elec, mgwe.twf, mgwe.ham, mgwe.particle_pool.getPool(), comm),
                         std::move(dmcdriver_input), *qmc_system,
                         MCPopulation(comm->size(), comm->rank(), &mgwe.pset_elec, &mgwe.twf, &mgwe.ham),
                         random_child_refs, comm);

  dmc_batched.process(node);
  dmc_batched.setStatus(test_project.currentMainRoot(), "", false);

  using Dmcbta     = testing::DMCBatchedTestAccessor;
  auto run_context = Dmcbta::startRun(dmc_batched);

  int global_step                                     = 0;
  int num_blocks                                      = qmcdriver_input.get_max_blocks();
  run_context.dmc_state.recalculate_properties_period = 100;
  run_context.dmc_state.is_recomputing_block          = true;
  auto& crowds                                        = Dmcbta::getCrowds(dmc_batched);
  app_log() << "crowds: " << crowds.size() << '\n';
  auto steps_per_block = Dmcbta::getStepsPerBlock(dmc_batched);
  for (UPtr<Crowd>& crowd : crowds)
    crowd->startBlock(steps_per_block);

  app_log() << "Beginning steps in test_DMCBatched_integration.cpp\n";
  app_log() << "For " << num_blocks << " Blocks.\n";

  int block = 0;
  while (block < num_blocks)
  {
    auto num_crowds = QMCDriverNew::determineNumCrowds(qmcdriver_input.get_num_crowds(), random_child_refs.size());

    for (int step = 0; step < steps_per_block; ++step, ++global_step)
    {
      run_context.dmc_state.step        = step;
      run_context.dmc_state.global_step = global_step;
      auto& context_for_steps           = Dmcbta::getContextForSteps(dmc_batched);

      //const int iter = block * steps_per_block + step;
      // walker_controller_->branch(iter, population_, iter == 0);
      // branch_engine_->updateParamAfterPopControl(walker_controller_->get_ensemble_property(),
      //                                            population_.get_golden_electrons().getTotalNum());
      // walker_controller_->setTrialEnergy(branch_engine_->getEtrial());
      auto& driver_timers = Dmcbta::getTimers(dmc_batched);
      for (int i_crowd = 0; i_crowd < num_crowds; ++i_crowd)
        Dmcbta::runDMCStep(dmc_batched, i_crowd, run_context.dmc_state, driver_timers,
                           Dmcbta::getDMCTimers(dmc_batched), context_for_steps, crowds);
    }
    Dmcbta::endBlock(dmc_batched, block);
    ++block;
  }
  Dmcbta::endRun(dmc_batched, block);
  outputManager.resume();

  // This gets a prefix based on the platform, the random numbers are
  // used in a sequence that is not "deterministic" between the
  // different platforms. So unless that can be fixed we need this
  std::string prefix{PLATFORM_PREFIX};
  // would be better to dig this out of the input but for expediency
  std::string log_name("rank_" + std::to_string(comm->rank()) + "_dmc_vem_test.dat");
  std::string log_dest_name(prefix + "_rank_" + std::to_string(comm->rank()) + "of" + std::to_string(comm->size()) +
                            "_dmc_vem_test.dat");

  if constexpr (generate_test_data)
  {
    copyRefFile(log_name, log_dest_name);
  }
  else
  {
    std::ifstream tf_stream(log_name.c_str());
    std::vector<std::string> tf_lines;
    for (std::string line; getline(tf_stream, line);)
      tf_lines.push_back(line);
    std::ifstream ref_stream(log_dest_name.c_str());
    std::vector<std::string> ref_lines;
    for (std::string line; getline(ref_stream, line);)
      ref_lines.push_back(line);
    CHECK(tf_lines == ref_lines);
  }
}

} // namespace qmcplusplus
