//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: DMC.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <functional>
#include <cassert>
#include <cmath>

#include "DMCBatched.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"
#include "Concurrency/ParallelExecutor.hpp"
#include "Concurrency/Info.hpp"
#include "Message/UniformCommunicateError.h"
#include "Message/CommOperators.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/RunTimeManager.h"
#include "Utilities/ProgressReportEngine.h"
#include "QMCDrivers/DMC/WalkerControl.h"
#include "QMCDrivers/SFNBranch.h"
#include "EstimatorInputDelegates.h"
#include "MemoryUsage.h"
#include "QMCWaveFunctions/TWFGrads.hpp"
#include "TauParams.hpp"

namespace qmcplusplus
{
using std::placeholders::_1;
using WP       = WalkerProperties::Indexes;
using PsiValue = TrialWaveFunction::PsiValue;

/** Constructor maintains proper ownership of input parameters
 *
 *  Note you must call the Base constructor before the derived class sets QMCType
 */
DMCBatched::DMCBatched(const ProjectData& project_data,
                       QMCDriverInput&& qmcdriver_input,
                       const std::optional<EstimatorManagerInput>& global_emi,
                       DMCDriverInput&& input,
                       WalkerConfigurations& wc,
                       MCPopulation&& pop,
                       Communicate* comm)
    : QMCDriverNew(project_data,
                   std::move(qmcdriver_input),
                   global_emi,
                   wc,
                   std::move(pop),
                   "DMCBatched::",
                   comm,
                   "DMCBatched",
                   std::bind(&DMCBatched::setNonLocalMoveHandler, this, _1)),
      dmcdriver_input_(input),
      dmc_timers_("DMCBatched::")
{}

DMCBatched::~DMCBatched() = default;

void DMCBatched::setNonLocalMoveHandler(QMCHamiltonian& golden_hamiltonian)
{
  golden_hamiltonian.setNonLocalMoves(dmcdriver_input_.get_non_local_move(), qmcdriver_input_.get_tau(),
                                      dmcdriver_input_.get_alpha(), dmcdriver_input_.get_gamma());
}

template<CoordsType CT>
void DMCBatched::advanceWalkers(const StateForThread& sft,
                                Crowd& crowd,
                                DriverTimers& timers,
                                DMCTimers& dmc_timers,
                                ContextForSteps& step_context,
                                bool recompute,
                                bool accumulate_this_step)
{
  auto& ps_dispatcher  = crowd.dispatchers_.ps_dispatcher_;
  auto& twf_dispatcher = crowd.dispatchers_.twf_dispatcher_;
  auto& ham_dispatcher = crowd.dispatchers_.ham_dispatcher_;

  auto& walkers = crowd.get_walkers();
  const RefVectorWithLeader<ParticleSet> walker_elecs(crowd.get_walker_elecs()[0], crowd.get_walker_elecs());
  const RefVectorWithLeader<TrialWaveFunction> walker_twfs(crowd.get_walker_twfs()[0], crowd.get_walker_twfs());
  const RefVectorWithLeader<QMCHamiltonian> walker_hamiltonians(crowd.get_walker_hamiltonians()[0],
                                                                crowd.get_walker_hamiltonians());

  timers.resource_timer.start();
  ResourceCollectionTeamLock<ParticleSet> pset_res_lock(crowd.getSharedResource().pset_res, walker_elecs);
  ResourceCollectionTeamLock<TrialWaveFunction> twfs_res_lock(crowd.getSharedResource().twf_res, walker_twfs);
  ResourceCollectionTeamLock<QMCHamiltonian> hams_res_lock(crowd.getSharedResource().ham_res, walker_hamiltonians);
  timers.resource_timer.stop();

  {
    ScopedTimer recompute_timer(dmc_timers.step_begin_recompute_timer);
    std::vector<bool> recompute_mask;
    recompute_mask.reserve(walkers.size());
    for (MCPWalker& awalker : walkers)
      if (awalker.wasTouched)
      {
        recompute_mask.push_back(true);
        awalker.wasTouched = false;
      }
      else
        recompute_mask.push_back(false);
    ps_dispatcher.flex_loadWalker(walker_elecs, walkers, recompute_mask, true);
    twf_dispatcher.flex_recompute(walker_twfs, walker_elecs, recompute_mask);
  }

  const int num_walkers   = crowd.size();
  auto& pset_leader       = walker_elecs.getLeader();
  const int num_particles = pset_leader.getTotalNum();

  MCCoords<CT> drifts(num_walkers), drifts_reverse(num_walkers);
  MCCoords<CT> walker_deltas(num_walkers * num_particles), deltas(num_walkers);
  TWFGrads<CT> grads_now(num_walkers), grads_new(num_walkers);

  //This generates an entire steps worth of deltas.
  makeGaussRandomWithEngine(walker_deltas, step_context.get_random_gen());

  std::vector<PsiValue> ratios(num_walkers, PsiValue(0.0));
  std::vector<RealType> log_gf(num_walkers, 0.0);
  std::vector<RealType> log_gb(num_walkers, 0.0);
  std::vector<RealType> prob(num_walkers, 0.0);

  // local list to handle accept/reject
  std::vector<bool> isAccepted;
  isAccepted.reserve(num_walkers);

  //save the old energies for branching needs.
  std::vector<FullPrecRealType> old_energies(num_walkers);
  for (int iw = 0; iw < num_walkers; ++iw)
    old_energies[iw] = walkers[iw].get().Properties(WP::LOCALENERGY);

  std::vector<RealType> rr_proposed(num_walkers, 0.0);
  std::vector<RealType> rr_accepted(num_walkers, 0.0);

  {
    ScopedTimer pbyp_local_timer(timers.movepbyp_timer);
    for (int ig = 0; ig < pset_leader.groups(); ++ig)
    {
      TauParams<RealType, CT> taus(sft.qmcdrv_input.get_tau(), sft.population.get_ptclgrp_inv_mass()[ig],
                                   sft.qmcdrv_input.get_spin_mass());

      twf_dispatcher.flex_prepareGroup(walker_twfs, walker_elecs, ig);

      for (int iat = pset_leader.first(ig); iat < pset_leader.last(ig); ++iat)
      {
        //This is very useful thing to be able to look at in the debugger
#ifndef NDEBUG
        std::vector<int> walkers_who_have_been_on_wire(num_walkers, 0);
        for (int iw = 0; iw < walkers.size(); ++iw)
        {
          walkers[iw].get().get_has_been_on_wire() ? walkers_who_have_been_on_wire[iw] = 1
                                                   : walkers_who_have_been_on_wire[iw] = 0;
        }
#endif
        //get deltas for this particle for all walkers
        walker_deltas.getSubset(iat * num_walkers, num_walkers, deltas);

        // only DMC does this
        // TODO: rr needs a real name
        std::vector<RealType> rr(num_walkers, 0.0);
        assert(rr.size() == deltas.positions.size());
        std::transform(deltas.positions.begin(), deltas.positions.end(), rr.begin(),
                       [t = taus.tauovermass](auto& delta_r) { return t * dot(delta_r, delta_r); });

        twf_dispatcher.flex_evalGrad(walker_twfs, walker_elecs, iat, grads_now);
        sft.drift_modifier.getDrifts(taus, grads_now, drifts);

        scaleBySqrtTau(taus, deltas);
        drifts += deltas;

// in DMC this was done here, changed to match VMCBatched pending factoring to common source
// if (rr > m_r2max)
// {
//   ++nRejectTemp;
//   continue;
// }
#ifndef NDEBUG
        for (int i = 0; i < rr.size(); ++i)
          assert(std::isfinite(rr[i]));
#endif

        ps_dispatcher.flex_makeMove(walker_elecs, iat, drifts);

        twf_dispatcher.flex_calcRatioGrad(walker_twfs, walker_elecs, iat, ratios, grads_new);

        computeLogGreensFunction(deltas, taus, log_gf);

        sft.drift_modifier.getDrifts(taus, grads_new, drifts_reverse);

        drifts_reverse += drifts;

        computeLogGreensFunction(drifts_reverse, taus, log_gb);

        auto checkPhaseChanged = [&sft](const PsiValue& ratio, int& is_reject) {
          if (ratio == PsiValue(0) || sft.branch_engine.phaseChanged(std::arg(ratio)))
            is_reject = 1;
          else
            is_reject = 0;
        };

        // Hopefully a phase change doesn't make any of these transformations fail.
        std::vector<int> rejects(num_walkers); // instead of std::vector<bool>
        for (int iw = 0; iw < num_walkers; ++iw)
        {
          checkPhaseChanged(ratios[iw], rejects[iw]);
          //This is just convenient to do here
          rr_proposed[iw] += rr[iw];
        }

        for (int iw = 0; iw < num_walkers; ++iw)
          prob[iw] = std::norm(ratios[iw]) * std::exp(log_gb[iw] - log_gf[iw]);

        isAccepted.clear();

        for (int iw = 0; iw < num_walkers; ++iw)
          if ((!rejects[iw]) && prob[iw] >= std::numeric_limits<RealType>::epsilon() &&
              step_context.get_random_gen()() < prob[iw])
          {
            crowd.incAccept();
            isAccepted.push_back(true);
            rr_accepted[iw] += rr[iw];
          }
          else
          {
            crowd.incReject();
            isAccepted.push_back(false);
          }

        twf_dispatcher.flex_accept_rejectMove(walker_twfs, walker_elecs, iat, isAccepted, true);

        ps_dispatcher.flex_accept_rejectMove<CT>(walker_elecs, iat, isAccepted);
      }
    }

    twf_dispatcher.flex_completeUpdates(walker_twfs);
    ps_dispatcher.flex_donePbyP(walker_elecs);
  }

  { // collect GL for KE.
    ScopedTimer buffer_local(timers.buffer_timer);
    twf_dispatcher.flex_evaluateGL(walker_twfs, walker_elecs, recompute);
    if (sft.qmcdrv_input.get_debug_checks() & DriverDebugChecks::CHECKGL_AFTER_MOVES)
      checkLogAndGL(crowd, "checkGL_after_moves");
    ps_dispatcher.flex_saveWalker(walker_elecs, walkers);
  }

  { // hamiltonian
    ScopedTimer ham_local(timers.hamiltonian_timer);

    std::vector<QMCHamiltonian::FullPrecRealType> new_energies(
        ham_dispatcher.flex_evaluateWithToperator(walker_hamiltonians, walker_twfs, walker_elecs));

    auto resetSigNLocalEnergy = [](MCPWalker& walker, TrialWaveFunction& twf, auto local_energy, auto rr_acc,
                                   auto rr_prop) {
      walker.resetProperty(twf.getLogPsi(), twf.getPhase(), local_energy, rr_acc, rr_prop, 1.0);
    };

    for (int iw = 0; iw < walkers.size(); ++iw)
    {
      resetSigNLocalEnergy(walkers[iw], walker_twfs[iw], new_energies[iw], rr_accepted[iw], rr_proposed[iw]);
      FullPrecRealType branch_weight = sft.branch_engine.branchWeight(new_energies[iw], old_energies[iw]);
      walkers[iw].get().Weight *= branch_weight;
      if (rr_proposed[iw] > 0)
        walkers[iw].get().Age = 0;
      else
        walkers[iw].get().Age++;
    }
  }

  { // estimator collectables
    ScopedTimer collectable_local(timers.collectables_timer);

    // evaluate non-physical hamiltonian elements
    for (int iw = 0; iw < walkers.size(); ++iw)
      walker_hamiltonians[iw].auxHevaluate(walker_elecs[iw], walkers[iw]);

    // save properties into walker
    for (int iw = 0; iw < walkers.size(); ++iw)
      walker_hamiltonians[iw].saveProperty(walkers[iw].get().getPropertyBase());
  }

  if (accumulate_this_step)
  {
    ScopedTimer est_timer(timers.estimators_timer);
    crowd.accumulate(step_context.get_random_gen());
  }

  { // T-moves
    ScopedTimer tmove_timer(dmc_timers.tmove_timer);

    const auto num_walkers = walkers.size();
    std::vector<int> walker_non_local_moves_accepted(num_walkers, 0);
    RefVector<MCPWalker> moved_nonlocal_walkers;
    RefVectorWithLeader<ParticleSet> moved_nonlocal_walker_elecs(crowd.get_walker_elecs()[0]);
    RefVectorWithLeader<TrialWaveFunction> moved_nonlocal_walker_twfs(crowd.get_walker_twfs()[0]);
    moved_nonlocal_walkers.reserve(num_walkers);
    moved_nonlocal_walker_elecs.reserve(num_walkers);
    moved_nonlocal_walker_twfs.reserve(num_walkers);

    for (int iw = 0; iw < walkers.size(); ++iw)
    {
      walker_non_local_moves_accepted[iw] = walker_hamiltonians[iw].makeNonLocalMoves(walker_elecs[iw]);

      if (walker_non_local_moves_accepted[iw] > 0)
      {
        crowd.incNonlocalAccept(walker_non_local_moves_accepted[iw]);
        moved_nonlocal_walkers.push_back(walkers[iw]);
        moved_nonlocal_walker_elecs.push_back(walker_elecs[iw]);
        moved_nonlocal_walker_twfs.push_back(walker_twfs[iw]);
      }
    }

    if (moved_nonlocal_walkers.size())
    {
      twf_dispatcher.flex_evaluateGL(moved_nonlocal_walker_twfs, moved_nonlocal_walker_elecs, false);
      if (sft.qmcdrv_input.get_debug_checks() & DriverDebugChecks::CHECKGL_AFTER_TMOVE)
        checkLogAndGL(crowd, "checkGL_after_tmove");
      ps_dispatcher.flex_saveWalker(moved_nonlocal_walker_elecs, moved_nonlocal_walkers);
    }
  }
}

template void DMCBatched::advanceWalkers<CoordsType::POS>(const StateForThread& sft,
                                                          Crowd& crowd,
                                                          DriverTimers& timers,
                                                          DMCTimers& dmc_timers,
                                                          ContextForSteps& step_context,
                                                          bool recompute,
                                                          bool accumulate_this_step);

template void DMCBatched::advanceWalkers<CoordsType::POS_SPIN>(const StateForThread& sft,
                                                               Crowd& crowd,
                                                               DriverTimers& timers,
                                                               DMCTimers& dmc_timers,
                                                               ContextForSteps& step_context,
                                                               bool recompute,
                                                               bool accumulate_this_step);

void DMCBatched::runDMCStep(int crowd_id,
                            const StateForThread& sft,
                            DriverTimers& timers,
                            DMCTimers& dmc_timers,
                            UPtrVector<ContextForSteps>& context_for_steps,
                            UPtrVector<Crowd>& crowds)
{
  Crowd& crowd = *(crowds[crowd_id]);

  if (crowd.size() == 0)
    return;

  auto& rng = context_for_steps[crowd_id]->get_random_gen();
  crowd.setRNGForHamiltonian(rng);

  const int max_steps  = sft.qmcdrv_input.get_max_steps();
  const IndexType step = sft.step;
  // Are we entering the the last step of a block to recompute at?
  const bool recompute_this_step  = (sft.is_recomputing_block && (step + 1) == max_steps);
  const bool accumulate_this_step = true;
  const bool spin_move            = sft.population.get_golden_electrons().isSpinor();
  if (spin_move)
    advanceWalkers<CoordsType::POS_SPIN>(sft, crowd, timers, dmc_timers, *context_for_steps[crowd_id],
                                         recompute_this_step, accumulate_this_step);
  else
    advanceWalkers<CoordsType::POS>(sft, crowd, timers, dmc_timers, *context_for_steps[crowd_id], recompute_this_step,
                                    accumulate_this_step);
}

void DMCBatched::process(xmlNodePtr node)
{
  print_mem("DMCBatched before initialization", app_log());
  try
  {
    QMCDriverNew::AdjustedWalkerCounts awc =
        adjustGlobalWalkerCount(*myComm, walker_configs_ref_.getActiveWalkers(), qmcdriver_input_.get_total_walkers(),
                                qmcdriver_input_.get_walkers_per_rank(), dmcdriver_input_.get_reserve(),
                                qmcdriver_input_.get_num_crowds());

    Base::initializeQMC(awc);
  }
  catch (const UniformCommunicateError& ue)
  {
    myComm->barrier_and_abort(ue.what());
  }

  {
    ReportEngine PRE("DMC", "resetUpdateEngines");
    Timer init_timer;
    // Here DMC loads "Ensemble of cloned MCWalkerConfigurations"
    // I'd like to do away with this method in DMCBatched.

    app_log() << "  Creating the branching engine and walker controler" << std::endl;
    const auto refE_update_scheme = dmcdriver_input_.get_refenergy_update_scheme();
    app_log() << "    Reference energy is updated using the "
              << (refE_update_scheme == DMCRefEnergyScheme::UNLIMITED_HISTORY ? "unlimited_history" : "limited_history")
              << " scheme" << std::endl;
    branch_engine_ =
        std::make_unique<SFNBranch>(qmcdriver_input_.get_tau(), dmcdriver_input_.get_feedback(), refE_update_scheme);
    branch_engine_->put(node);

    walker_controller_ = std::make_unique<WalkerControl>(myComm, Random, dmcdriver_input_.get_reconfiguration());
    walker_controller_->setMinMax(population_.get_num_global_walkers(), 0);
    walker_controller_->start();
    walker_controller_->put(node);

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
}

bool DMCBatched::run()
{
  IndexType num_blocks = qmcdriver_input_.get_max_blocks();

  estimator_manager_->startDriverRun();
  StateForThread dmc_state(qmcdriver_input_, dmcdriver_input_, *drift_modifier_, *branch_engine_, population_);

  LoopTimer<> dmc_loop;
  RunTimeControl<> runtimeControl(run_time_manager, project_data_.getMaxCPUSeconds(), project_data_.getTitle(),
                                  myComm->rank() == 0);

  { // walker initialization
    ScopedTimer local_timer(timers_.init_walkers_timer);
    ParallelExecutor<> section_start_task;
    section_start_task(crowds_.size(), initialLogEvaluation, std::ref(crowds_), std::ref(step_contexts_));

    FullPrecRealType energy, variance;
    population_.measureGlobalEnergyVariance(*myComm, energy, variance);
    // false indicates we do not support kill at node crossings.
    branch_engine_->initParam(population_, energy, variance, dmcdriver_input_.get_reconfiguration(), false);
    walker_controller_->setTrialEnergy(branch_engine_->getEtrial());

    print_mem("DMCBatched after initialLogEvaluation", app_summary());
    if (qmcdriver_input_.get_measure_imbalance())
      measureImbalance("InitialLogEvaluation");
  }

  // this barrier fences all previous load imbalance. Avoid block 0 timing pollution.
  myComm->barrier();

  ScopedTimer local_timer(timers_.production_timer);
  ParallelExecutor<> crowd_task;

  for (int block = 0; block < num_blocks; ++block)
  {
    {
      ScopeGuard<LoopTimer<>> dmc_local_timer(dmc_loop);
      estimator_manager_->startBlock(qmcdriver_input_.get_max_steps());

      dmc_state.recalculate_properties_period = (qmc_driver_mode_[QMC_UPDATE_MODE])
          ? qmcdriver_input_.get_recalculate_properties_period()
          : (qmcdriver_input_.get_max_blocks() + 1) * qmcdriver_input_.get_max_steps();
      dmc_state.is_recomputing_block          = qmcdriver_input_.get_blocks_between_recompute()
                   ? (1 + block) % qmcdriver_input_.get_blocks_between_recompute() == 0
                   : false;

      for (UPtr<Crowd>& crowd : crowds_)
        crowd->startBlock(qmcdriver_input_.get_max_steps());

      for (int step = 0; step < qmcdriver_input_.get_max_steps(); ++step)
      {
        ScopedTimer local_timer(timers_.run_steps_timer);
        dmc_state.step = step;
        crowd_task(crowds_.size(), runDMCStep, dmc_state, timers_, dmc_timers_, std::ref(step_contexts_),
                   std::ref(crowds_));

        {
          const int iter = block * qmcdriver_input_.get_max_steps() + step;
          walker_controller_->branch(iter, population_, iter == 0);
          branch_engine_->updateParamAfterPopControl(walker_controller_->get_ensemble_property(),
                                                     population_.get_golden_electrons().getTotalNum());
          walker_controller_->setTrialEnergy(branch_engine_->getEtrial());
        }

        population_.redistributeWalkers(crowds_);
      }
      print_mem("DMCBatched after a block", app_debug_stream());
      if (qmcdriver_input_.get_measure_imbalance())
        measureImbalance("Block " + std::to_string(block));
      endBlock();
      recordBlock(block);
    }

    bool stop_requested = false;
    // Rank 0 decides whether the time limit was reached
    if (!myComm->rank())
      stop_requested = runtimeControl.checkStop(dmc_loop);
    myComm->bcast(stop_requested);

    if (stop_requested)
    {
      if (!myComm->rank())
        app_log() << runtimeControl.generateStopMessage("DMCBatched", block);
      run_time_manager.markStop();
      break;
    }
  }

  branch_engine_->printStatus();

  print_mem("DMCBatched ends", app_log());

  estimator_manager_->stopDriverRun();

  return finalize(num_blocks, true);
}

} // namespace qmcplusplus
