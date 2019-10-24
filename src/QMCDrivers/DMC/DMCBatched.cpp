//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: DMC.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <functional>

#include "QMCDrivers/DMC/DMCBatched.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"
#include "Concurrency/TasksOneToOne.hpp"
#include "Concurrency/Info.hpp"
#include "Utilities/RunTimeManager.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
using std::placeholders::_1;

// clang-format off
/** Constructor maintains proper ownership of input parameters
   */
DMCBatched::DMCBatched(QMCDriverInput&& qmcdriver_input,
                       DMCDriverInput&& input,
                       MCPopulation& pop,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       WaveFunctionPool& wf_pool,
                       Communicate* comm)
    : QMCDriverNew(std::move(qmcdriver_input), pop, psi, h, wf_pool,
                   "DMCBatched::", comm,
                   std::bind(&DMCBatched::setNonLocalMoveHandler, this, _1)),
      dmcdriver_input_(input)
{
  QMCType = "DMCBatched";
}
// clang-format on

QMCTraits::IndexType DMCBatched::calc_default_local_walkers(IndexType walkers_per_rank)
{
  checkNumCrowdsLTNumThreads();
  int num_threads(Concurrency::maxThreads<>());
  IndexType rw = walkers_per_rank; //qmcdriver_input_.get_walkers_per_rank();
  if (num_crowds_ == 0)
    num_crowds_ = std::min(num_threads, rw);
  walkers_per_crowd_ = (rw % num_crowds_) ? rw / num_crowds_ + 1 : rw / num_crowds_;

  IndexType local_walkers = walkers_per_crowd_ * num_crowds_;
  population_.set_num_local_walkers(local_walkers);
  population_.set_num_global_walkers(local_walkers * population_.get_num_ranks());
  if (rw != qmcdriver_input_.get_walkers_per_rank())
    app_warning() << "DMCBatched driver has adjusted walkers per rank to: " << local_walkers << '\n';

  app_log() << "DMCBatched walkers per crowd " << walkers_per_crowd_ << std::endl;
  return local_walkers;
}

void DMCBatched::setNonLocalMoveHandler(QMCHamiltonian& golden_hamiltonian)
{
  golden_hamiltonian.setNonLocalMoves(dmcdriver_input_.get_non_local_move(), qmcdriver_input_.get_tau(),
                                      dmcdriver_input_.get_alpha(), dmcdriver_input_.get_gamma());
}

void DMCBatched::resetUpdateEngines()
{
  ReportEngine PRE("DMC", "resetUpdateEngines");
  Timer init_timer;
  // Here DMC loads "Ensemble of cloned MCWalkerConfigurations"
  int nw_multi = branch_engine_->initWalkerController(population_, dmcdriver_input_.get_reconfiguration(), false);
  RefVector<MCPWalker> walkers(convertUPtrToRefVector(population_.get_walkers()));

  branch_engine_->checkParameters(population_.get_num_global_walkers(), walkers);

  std::ostringstream o;
  if (dmcdriver_input_.get_reconfiguration())
    o << "  Fixed population using reconfiguration method\n";
  else
    o << "  Fluctuating population\n";
  o << "  Persistent walkers are killed after " << dmcdriver_input_.get_max_age() << " MC sweeps\n";
  o << "  BranchInterval = " << dmcdriver_input_.get_branch_interval() << "\n";
  o << "  Steps per block = " << qmcdriver_input_.get_max_steps() << "\n";
  o << "  Number of blocks = " << qmcdriver_input_.get_max_blocks() << "\n";
  app_log() << o.str() << std::endl;

  app_log() << "  DMC Engine Initialization = " << init_timer.elapsed() << " secs" << std::endl;
}

void DMCBatched::advanceWalkers(const StateForThread& sft,
                                Crowd& crowd,
                                DriverTimers& timers,
                                //                                DMCTimers& dmc_timers,
                                ContextForSteps& step_context,
                                bool recompute)
{
  timers.buffer_timer.start();
  crowd.loadWalkers();

  constexpr RealType mhalf(-0.5);

  int naccept(0);
  int nreject(0);
  int nnode_crossing(0);

  // Consider favoring lambda followed by for over walkers
  // more compact, descriptive and less error prone.
  auto& walker_twfs      = crowd.get_walker_twfs();
  auto& walkers          = crowd.get_walkers();
  auto& walker_elecs     = crowd.get_walker_elecs();
  auto copyTWFFromBuffer = [](TrialWaveFunction& twf, ParticleSet& pset, MCPWalker& walker) {
    twf.copyFromBuffer(pset, walker.DataSet);
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    copyTWFFromBuffer(walker_twfs[iw], walker_elecs[iw], walkers[iw]);
  timers.buffer_timer.stop();

  timers.movepbyp_timer.start();
  int num_walkers = crowd.size();
  //This generates an entire steps worth of deltas.
  step_context.nextDeltaRs();
  auto it_delta_r = step_context.deltaRsBegin();
  std::vector<PosType> drifts(num_walkers);

  // local list to handle accept/reject
  std::vector<std::reference_wrapper<ParticleSet>> elec_accept_list, elec_reject_list;
  std::vector<std::reference_wrapper<TrialWaveFunction>> twf_accept_list, twf_reject_list;
  elec_accept_list.reserve(num_walkers);
  elec_reject_list.reserve(num_walkers);
  twf_accept_list.reserve(num_walkers);
  twf_reject_list.reserve(num_walkers);

  //copy the old energy
  std::vector<FullPrecRealType> old_walker_energies(num_walkers);
  auto setOldEnergies = [](MCPWalker& walker, FullPrecRealType& old_walker_energy) {
    old_walker_energy = walker.Properties(LOCALENERGY);
  };
  for (int iw = 0; iw < num_walkers; ++iw)
    setOldEnergies(walkers[iw], old_walker_energies[iw]);
  std::vector<FullPrecRealType> new_walker_energies(old_walker_energies);

  std::vector<RealType> gf_acc(num_walkers, 1.0);
  std::vector<RealType> rr_proposed(num_walkers, 0.0);
  std::vector<RealType> rr_accepted(num_walkers, 0.0);

  //  dmc_timers_.dmc_movePbyP.start();

  std::vector<int> did_walker_move(num_walkers, 0);

  for (int ig = 0; ig < step_context.get_num_groups(); ++ig)
  {
    RealType tauovermass = sft.qmcdrv_input.get_tau() * sft.population.get_ptclgrp_inv_mass()[ig];
    RealType oneover2tau = 0.5 / (tauovermass);
    RealType sqrttau     = std::sqrt(tauovermass);
    int start_index      = step_context.getPtclGroupStart(ig);
    int end_index        = step_context.getPtclGroupEnd(ig);
    for (int iat = start_index; iat < end_index; ++iat)
    {
      // Correct for DMC as well as VMC?
      crowd.clearResults();
      ParticleSet::flex_setActive(crowd.get_walker_elecs(), iat);
      auto delta_r_start = it_delta_r + iat * num_walkers;
      auto delta_r_end   = delta_r_start + num_walkers;

      //get the displacement
      TrialWaveFunction::flex_evalGrad(crowd.get_walker_twfs(), crowd.get_walker_elecs(), iat, crowd.get_grads_now());
      sft.drift_modifier.getDrifts(tauovermass, crowd.get_grads_now(), drifts);

      std::transform(drifts.begin(), drifts.end(), delta_r_start, drifts.begin(),
                     [sqrttau](PosType& drift, PosType& delta_r) { return drift + (sqrttau * delta_r); });

      // only DMC does this
      // TODO: rr needs a real name
      std::vector<RealType> rr(num_walkers);
      std::transform(delta_r_start, delta_r_end, rr.begin(),
                     [tauovermass](auto& delta_r) { return tauovermass * dot(delta_r, delta_r); });

      // in DMC this was done here, changed to match VMCBatched pending factoring to common source
      // if (rr > m_r2max)
      // {
      //   ++nRejectTemp;
      //   continue;
      // }
      auto elecs = crowd.get_walker_elecs();
      ParticleSet::flex_makeMove(crowd.get_walker_elecs(), iat, drifts);

      TrialWaveFunction::flex_ratioGrad(crowd.get_walker_twfs(), crowd.get_walker_elecs(), iat, crowd.get_ratios(),
                                        crowd.get_grads_new());

      // This lambda is not nested thread safe due to the nreject, nnode_crossing updates
      auto checkPhaseChanged = [&sft, &iat, &crowd, &nnode_crossing](TrialWaveFunction& twf, ParticleSet& elec,
                                                                     int& is_reject) {
        if (sft.branch_engine.phaseChanged(twf.getPhaseDiff()))
        {
          crowd.incReject();
          ++nnode_crossing;
          elec.rejectMove(iat);
          twf.rejectMove(iat);
          is_reject = 1;
        }
        else
          is_reject = 0;
      };

      // Hopefully a phase change doesn't make any of these transformations fail.

      std::vector<int> rejects(num_walkers);
      for (int iw = 0; iw < num_walkers; ++iw)
      {
        checkPhaseChanged(walker_twfs[iw], walker_elecs[iw], rejects[iw]);
        //This is just convenient to do here
        rr_proposed[iw] += rr[iw];
      }
      std::transform(delta_r_start, delta_r_end, crowd.get_log_gf().begin(),
                     [mhalf](auto& delta_r) { return mhalf * dot(delta_r, delta_r); });

      sft.drift_modifier.getDrifts(tauovermass, crowd.get_grads_new(), drifts);

      std::transform(crowd.beginElectrons(), crowd.endElectrons(), drifts.begin(), drifts.begin(),
                     [iat](auto& elecs, auto& drift) { return elecs.get().R[iat] - elecs.get().activePos - drift; });

      std::transform(drifts.begin(), drifts.end(), crowd.get_log_gb().begin(),
                     [oneover2tau](auto& drift) { return -oneover2tau * dot(drift, drift); });

      std::transform(crowd.get_ratios().begin(), crowd.get_ratios().end(), crowd.get_prob().begin(),
                     [](auto ratio) { return std::norm(ratio); });

      twf_accept_list.clear();
      twf_reject_list.clear();
      elec_accept_list.clear();
      elec_reject_list.clear();


      for (int iw = 0; iw < num_walkers; ++iw)
      {
        auto prob   = crowd.get_prob()[iw];
        auto log_gf = crowd.get_log_gf()[iw];
        auto log_gb = crowd.get_log_gb()[iw];

        if ((!rejects[iw]) && prob >= std::numeric_limits<RealType>::epsilon() &&
            step_context.get_random_gen()() < prob * std::exp(log_gb - log_gf))
        {
          did_walker_move[iw] += 1;
          crowd.incAccept();
          twf_accept_list.push_back(crowd.get_walker_twfs()[iw]);
          elec_accept_list.push_back(crowd.get_walker_elecs()[iw]);
          rr_accepted[iw] += rr[iw];
          gf_acc[iw] *= prob;
        }
        else
        {
          crowd.incReject();
          twf_reject_list.push_back(crowd.get_walker_twfs()[iw]);
          elec_reject_list.push_back(crowd.get_walker_elecs()[iw]);
        }
      }

      TrialWaveFunction::flex_acceptMove(twf_accept_list, elec_accept_list, iat);
      TrialWaveFunction::flex_rejectMove(twf_reject_list, iat);

      ParticleSet::flex_acceptMove(elec_accept_list, iat);
      ParticleSet::flex_rejectMove(elec_reject_list, iat);
    }
  }
  std::for_each(crowd.get_walker_twfs().begin(), crowd.get_walker_twfs().end(),
                [](auto& twf) { twf.get().completeUpdates(); });

  ParticleSet::flex_donePbyP(crowd.get_walker_elecs());
  //dmc_timers.dmc_movePbyP.stop();
  timers.movepbyp_timer.stop();

  //To use the flex interfaces we have to build RefVectors for walker that moved and walkers that didn't

  auto& walker_hamiltonians  = crowd.get_walker_hamiltonians();
  auto& walker_mcp_wfbuffers = crowd.get_mcp_wfbuffers();

  int num_moved = 0;
  for (int iw = 0; iw < num_walkers; ++iw)
  {
    if (did_walker_move[iw] > 0)
      num_moved++;
  }

  DMCPerWalkerRefs moved(num_moved);
  DMCPerWalkerRefs stalled(num_walkers - num_moved);

  for (int iw = 0; iw < num_walkers; ++iw)
  {
    if (did_walker_move[iw] > 0)
    {
      moved.walkers.push_back(walkers[iw]);
      moved.walker_twfs.push_back(walker_twfs[iw]);
      moved.walker_hamiltonians.push_back(walker_hamiltonians[iw]);
      moved.walker_elecs.push_back(walker_elecs[iw]);
      moved.walker_mcp_wfbuffers.push_back(walker_mcp_wfbuffers[iw]);
      moved.old_energies.push_back(old_walker_energies[iw]);
      moved.new_energies.push_back(new_walker_energies[iw]);
    }
    else
    {
      stalled.walkers.push_back(walkers[iw]);
      stalled.walker_twfs.push_back(walker_twfs[iw]);
      stalled.walker_hamiltonians.push_back(walker_hamiltonians[iw]);
      stalled.walker_elecs.push_back(walker_elecs[iw]);
      stalled.walker_mcp_wfbuffers.push_back(walker_mcp_wfbuffers[iw]);
      stalled.old_energies.push_back(old_walker_energies[iw]);
      stalled.new_energies.push_back(new_walker_energies[iw]);
      stalled.gf_accs.push_back(gf_acc[iw]);
    }
  }

  if (moved.walkers.size() > 0)
  {
    timers.buffer_timer.start();
    TrialWaveFunction::flex_updateBuffer(moved.walker_twfs, moved.walker_elecs, moved.walker_mcp_wfbuffers);

    ParticleSet::flex_saveWalker(moved.walker_elecs, moved.walkers);
    timers.buffer_timer.stop();
    timers.hamiltonian_timer.start();
    // std::vector<QMCHamiltonian::FullPrecRealType> local_energies(QMCHamiltonian::flex_evaluate(moved.walker_hamiltonians, moved.walker_elecs));
    std::vector<QMCHamiltonian::FullPrecRealType> local_energies(
        QMCHamiltonian::flex_evaluateWithToperator(moved.walker_hamiltonians, moved.walker_elecs));
    timers.hamiltonian_timer.stop();

    auto resetSigNLocalEnergy = [](MCPWalker& walker, TrialWaveFunction& twf, auto& local_energy) {
      walker.resetProperty(twf.getLogPsi(), twf.getPhase(), local_energy);
    };
    for (int iw = 0; iw < moved.walkers.size(); ++iw)
    {
      resetSigNLocalEnergy(moved.walkers[iw], moved.walker_twfs[iw], local_energies[iw]);
      moved.walkers[iw].get().Weight *=
          sft.branch_engine.branchWeightBare(moved.new_energies[iw], moved.old_energies[iw]);
    }

    timers.collectables_timer.start();
    auto evaluateNonPhysicalHamiltonianElements = [](QMCHamiltonian& ham, ParticleSet& pset, MCPWalker& walker) {
      ham.auxHevaluate(pset, walker);
    };
    for (int iw = 0; iw < moved.walkers.size(); ++iw)
      evaluateNonPhysicalHamiltonianElements(moved.walker_hamiltonians[iw], moved.walker_elecs[iw], moved.walkers[iw]);

    auto savePropertiesIntoWalker = [](QMCHamiltonian& ham, MCPWalker& walker) {
      ham.saveProperty(walker.getPropertyBase());
    };
    for (int iw = 0; iw < moved.walkers.size(); ++iw)
      savePropertiesIntoWalker(moved.walker_hamiltonians[iw], moved.walkers[iw]);
    timers.collectables_timer.stop();
  }
  for (int iw = 0; iw < stalled.walkers.size(); ++iw)
  {
    MCPWalker& stalled_walker = stalled.walkers[iw];
    stalled_walker.Age++;
    stalled_walker.Properties(R2ACCEPTED) = 0.0;
    RealType wtmp                         = stalled_walker.Weight;
    // TODO: fix this walker.Weight twiddle for rejectedMove
    stalled_walker.Weight                      = 0.0;
    QMCHamiltonian& stalled_walker_hamiltonian = stalled.walker_hamiltonians[iw];
    ParticleSet& stalled_particle_set          = stalled.walker_elecs[iw];
    stalled_walker_hamiltonian.rejectedMove(stalled_particle_set, stalled_walker);
    stalled_walker.Weight                       = wtmp;
    FullPrecRealType& stalled_new_walker_energy = stalled.new_energies[iw];
    FullPrecRealType& stalled_old_walker_energy = stalled.old_energies[iw];
    stalled_new_walker_energy                   = stalled_old_walker_energy;
    RealType& gf_acc                            = stalled.gf_accs[iw];
    gf_acc                                      = 1.0;
    stalled_walker.Weight *= sft.branch_engine.branchWeight(stalled_new_walker_energy, stalled_old_walker_energy);
  }

  QMCHamiltonian& db_hamiltonian = walker_hamiltonians[0].get();


  //myTimers[DMC_tmoves]->start();
  std::vector<int> walker_non_local_moves_accepted(
      QMCHamiltonian::flex_makeNonLocalMoves(walker_hamiltonians, walker_elecs));


  // could be premature optimization
  int num_moved_nonlocal   = 0;
  int total_moved_nonlocal = 0;
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    if (walker_non_local_moves_accepted[iw] > 0)
    {
      num_moved_nonlocal++;
      total_moved_nonlocal += walker_non_local_moves_accepted[iw];
    }
  }

  if (num_moved_nonlocal > 0)
  {
    DMCPerWalkerRefs moved_nonlocal(num_moved_nonlocal);

    for (int iw = 0; iw < moved.walkers.size(); ++iw)
    {
      if (walker_non_local_moves_accepted[iw] > 0)
      {
        moved_nonlocal.walkers.push_back(walkers[iw]);
        moved_nonlocal.walker_twfs.push_back(walker_twfs[iw]);
        moved_nonlocal.walker_elecs.push_back(walker_elecs[iw]);
        moved_nonlocal.walker_mcp_wfbuffers.push_back(walker_mcp_wfbuffers[iw]);
      }
    }
    TrialWaveFunction::flex_updateBuffer(moved_nonlocal.walker_twfs, moved_nonlocal.walker_elecs,
                                         moved_nonlocal.walker_mcp_wfbuffers);
    ParticleSet::flex_saveWalker(moved_nonlocal.walker_elecs, moved_nonlocal.walkers);
    crowd.incNonlocalAccept(total_moved_nonlocal);
  }
  setMultiplicities(sft.dmcdrv_input, walkers, step_context.get_random_gen());
}

void DMCBatched::setMultiplicities(const DMCDriverInput& dmcdriver_input,
                                   RefVector<MCPWalker>& walkers,
                                   RandomGenerator_t& rng)
{
  auto setMultiplicity = [&dmcdriver_input, &rng](MCPWalker& walker) {
    constexpr RealType onehalf(0.5);
    constexpr RealType cone(1);
    RealType M = walker.Weight;
    if (walker.Age > dmcdriver_input.get_max_age())
      M = std::min(onehalf, M);
    else if (walker.Age > 0)
      M = std::min(cone, M);
    walker.Multiplicity = M + rng();
  };
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    setMultiplicity(walkers[iw]);
  }
}


void DMCBatched::runDMCStep(int crowd_id,
                            const StateForThread& sft,
                            DriverTimers& timers,
                            //                            DMCTimers& dmc_timers,
                            UPtrVector<ContextForSteps>& context_for_steps,
                            UPtrVector<Crowd>& crowds)
{
  Crowd& crowd                         = *(crowds[crowd_id]);
  crowd.setRNGForHamiltonian(context_for_steps[crowd_id]->get_random_gen());

  int max_steps = sft.qmcdrv_input.get_max_steps();
  // This is migraine inducing here and in the original driver, I believe they are the same in
  // VMC(Batched)/DMC(Batched) needs another check and unit test
  bool is_recompute_block =
      sft.recomputing_blocks ? (1 + sft.block) % sft.qmcdrv_input.get_blocks_between_recompute() == 0 : false;
  IndexType step           = sft.step;
  bool recompute_this_step = (is_recompute_block && (step + 1) == max_steps);
  advanceWalkers(sft, crowd, timers, *context_for_steps[crowd_id], recompute_this_step);
}

bool DMCBatched::run()
{
  resetUpdateEngines();
  IndexType num_blocks = qmcdriver_input_.get_max_blocks();

  estimator_manager_->setCollectionMode(true);
  estimator_manager_->start(num_blocks);
  StateForThread dmc_state(qmcdriver_input_, dmcdriver_input_, *drift_modifier_, *branch_engine_, population_);

  LoopTimer dmc_loop;

  int sample = 0;
  RunTimeControl runtimeControl(RunTimeManager, MaxCPUSecs);

  { // walker initialization
    ScopedTimer local_timer(&(timers_.init_walkers_timer));
    TasksOneToOne<> section_start_task(num_crowds_);
    section_start_task(initialLogEvaluation, std::ref(crowds_), std::ref(step_contexts_));
  }


  TasksOneToOne<> crowd_task(num_crowds_);

  for (int block = 0; block < num_blocks; ++block)
  {
    dmc_loop.start();
    estimator_manager_->startBlock(qmcdriver_input_.get_max_steps());

    dmc_state.recalculate_properties_period = (qmc_driver_mode_[QMC_UPDATE_MODE])
        ? qmcdriver_input_.get_recalculate_properties_period()
        : (qmcdriver_input_.get_max_blocks() + 1) * qmcdriver_input_.get_max_steps();
    dmc_state.recomputing_blocks = qmcdriver_input_.get_blocks_between_recompute();

    for (auto& crowd : crowds_)
      crowd->startBlock(qmcdriver_input_.get_max_steps());

    for (int step = 0; step < qmcdriver_input_.get_max_steps(); ++step)
    {
      ScopedTimer local_timer(&(timers_.run_steps_timer));
      dmc_state.step = step;
      crowd_task(runDMCStep, dmc_state, timers_, std::ref(step_contexts_), std::ref(crowds_));

      branch_engine_->branch(step, crowds_, population_);
      for (auto& crowd_ptr : crowds_)
        crowd_ptr->clearWalkers();
      population_.distributeWalkers(crowds_.begin(), crowds_.end(), walkers_per_crowd_);
    }

    RefVector<ScalarEstimatorBase> all_scalar_estimators;
    FullPrecRealType total_block_weight = 0.0;
    FullPrecRealType total_accept_ratio = 0.0;
    // Collect all the ScalarEstimatorsFrom EMCrowds
    for (const UPtr<Crowd>& crowd : crowds_)
    {
      auto crowd_sc_est = crowd->get_estimator_manager_crowd().get_scalar_estimators();
      all_scalar_estimators.insert(all_scalar_estimators.end(), std::make_move_iterator(crowd_sc_est.begin()),
                                   std::make_move_iterator(crowd_sc_est.end()));
      total_block_weight += crowd->get_estimator_manager_crowd().get_block_weight();
      total_accept_ratio += crowd->get_accept_ratio();
    }
    // Should this be adjusted if crowds have different
    total_accept_ratio /= crowds_.size();
    estimator_manager_->collectScalarEstimators(all_scalar_estimators, population_.get_num_local_walkers(),
                                                total_block_weight);
    // TODO: should be accept rate for block
    estimator_manager_->stopBlockNew(total_accept_ratio);
  }
  return false;
}

} // namespace qmcplusplus
