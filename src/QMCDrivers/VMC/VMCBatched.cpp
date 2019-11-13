//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: VMC.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/VMC/VMCBatched.h"
#include "Concurrency/TasksOneToOne.hpp"
#include "Concurrency/Info.hpp"
#include "Utilities/RunTimeManager.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus
{
/** Constructor maintains proper ownership of input parameters
   */
VMCBatched::VMCBatched(QMCDriverInput&& qmcdriver_input,
                       VMCDriverInput&& input,
                       MCPopulation& pop,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       WaveFunctionPool& ppool,
                       Communicate* comm)
    : QMCDriverNew(std::move(qmcdriver_input), pop, psi, h, ppool, "VMCBatched::", comm), vmcdriver_input_(input)
{
  QMCType = "VMCBatched";
  // qmc_driver_mode.set(QMC_UPDATE_MODE, 1);
  // qmc_driver_mode.set(QMC_WARMUP, 0);
}

VMCBatched::IndexType VMCBatched::calc_default_local_walkers(IndexType walkers_per_rank)
{
  checkNumCrowdsLTNumThreads();
  int num_threads(Concurrency::maxThreads<>());
  if (num_crowds_ == 0)
    num_crowds_ = std::min(num_threads, walkers_per_rank);

  if (walkers_per_rank < num_crowds_)
    walkers_per_rank = num_crowds_;
  walkers_per_crowd_ =
      (walkers_per_rank % num_crowds_) ? walkers_per_rank / num_crowds_ + 1 : walkers_per_rank / num_crowds_;
  IndexType local_walkers = walkers_per_crowd_ * num_crowds_;
  population_.set_num_local_walkers(local_walkers);
  population_.set_num_global_walkers(local_walkers * population_.get_num_ranks());
  if (walkers_per_rank != qmcdriver_input_.get_walkers_per_rank())
    app_warning() << "VMCBatched driver has adjusted walkers per rank to: " << local_walkers << '\n';

  if (vmcdriver_input_.get_samples() >= 0 || vmcdriver_input_.get_samples_per_thread() >= 0 ||
      vmcdriver_input_.get_steps_between_samples() >= 0)
    app_warning() << "VMCBatched currently ignores samples and samplesperthread\n";

  if (local_walkers != walkers_per_rank)
    app_warning() << "VMCBatched changed the number of walkers to " << local_walkers << ". User input was "
                  << walkers_per_rank << std::endl;

  app_log() << "VMCBatched walkers per crowd " << walkers_per_crowd_ << std::endl;
  // TODO: Simplify samples, samples per thread etc in the unified driver
  // see logic in original VMC.cpp
  return local_walkers;
}

void VMCBatched::advanceWalkers(const StateForThread& sft,
                                Crowd& crowd,
                                QMCDriverNew::DriverTimers& timers,
                                ContextForSteps& step_context,
                                bool recompute)
{
  timers.buffer_timer.start();
  crowd.loadWalkers();

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
  const int num_walkers = crowd.size();
  // Note std::vector<bool> is not like the rest of stl.
  std::vector<bool> moved(num_walkers, false);
  constexpr RealType mhalf(-0.5);
  const bool use_drift = sft.vmcdrv_input.get_use_drift();
  std::vector<TrialWaveFunction::GradType> grads_now(num_walkers);
  std::vector<TrialWaveFunction::GradType> grads_new(num_walkers);
  std::vector<TrialWaveFunction::PsiValueType> ratios(num_walkers);

  std::vector<PosType> drifts(num_walkers);
  std::vector<RealType> log_gf(num_walkers);
  std::vector<RealType> log_gb(num_walkers);
  std::vector<RealType> prob(num_walkers);

  // local list to handle accept/reject
  std::vector<std::reference_wrapper<ParticleSet>> elec_accept_list, elec_reject_list;
  std::vector<std::reference_wrapper<TrialWaveFunction>> twf_accept_list, twf_reject_list;
  elec_accept_list.reserve(num_walkers);
  elec_reject_list.reserve(num_walkers);
  twf_accept_list.reserve(num_walkers);
  twf_reject_list.reserve(num_walkers);

  for (int sub_step = 0; sub_step < sft.qmcdrv_input.get_sub_steps(); sub_step++)
  {
    //This generates an entire steps worth of deltas.
    step_context.nextDeltaRs();

    // up and down electrons are "species" within qmpack
    for (int ig = 0; ig < step_context.get_num_groups(); ++ig) //loop over species
    {
      RealType tauovermass = sft.qmcdrv_input.get_tau() * sft.population.get_ptclgrp_inv_mass()[ig];
      RealType oneover2tau = 0.5 / (tauovermass);
      RealType sqrttau     = std::sqrt(tauovermass);
      int start_index      = step_context.getPtclGroupStart(ig);
      int end_index        = step_context.getPtclGroupEnd(ig);
      for (int iat = start_index; iat < end_index; ++iat)
      {
        ParticleSet::flex_setActive(crowd.get_walker_elecs(), iat);
        // step_context.deltaRsBegin returns an iterator to a flat series of PosTypes
        // fastest in walkers then particles
        auto delta_r_start = step_context.deltaRsBegin() + iat * num_walkers;
        auto delta_r_end   = delta_r_start + num_walkers;

        if (use_drift)
        {
          TrialWaveFunction::flex_evalGrad(crowd.get_walker_twfs(), crowd.get_walker_elecs(), iat, grads_now);
          sft.drift_modifier.getDrifts(tauovermass, grads_now, drifts);

          std::transform(drifts.begin(), drifts.end(), delta_r_start, drifts.begin(),
                         [sqrttau](const PosType& drift, const PosType& delta_r) {
                           return drift + (sqrttau * delta_r);
                         });
        }
        else
        {
          std::transform(delta_r_start, delta_r_end, drifts.begin(),
                         [sqrttau](const PosType& delta_r) { return sqrttau * delta_r; });
        }

        ParticleSet::flex_makeMove(crowd.get_walker_elecs(), iat, drifts);

        // This is inelegant
        if (use_drift)
        {
          TrialWaveFunction::flex_ratioGrad(crowd.get_walker_twfs(), crowd.get_walker_elecs(), iat, ratios, grads_new);
          std::transform(delta_r_start, delta_r_end, log_gf.begin(),
                         [mhalf](const PosType& delta_r) { return mhalf * dot(delta_r, delta_r); });

          sft.drift_modifier.getDrifts(tauovermass, grads_new, drifts);

          std::transform(crowd.beginElectrons(), crowd.endElectrons(), drifts.begin(), drifts.begin(),
                         [iat](const ParticleSet& elecs, const PosType& drift) {
                           return elecs.R[iat] - elecs.activePos - drift;
                         });

          std::transform(drifts.begin(), drifts.end(), log_gb.begin(),
                         [oneover2tau](const PosType& drift) { return -oneover2tau * dot(drift, drift); });
        }
        else
        {
          TrialWaveFunction::flex_calcRatio(crowd.get_walker_twfs(), crowd.get_walker_elecs(), iat, ratios);
        }

        std::transform(ratios.begin(), ratios.end(), prob.begin(), [](auto ratio) { return std::norm(ratio); });

        twf_accept_list.clear();
        twf_reject_list.clear();
        elec_accept_list.clear();
        elec_reject_list.clear();

        for (int i_accept = 0; i_accept < num_walkers; ++i_accept)
          if (prob[i_accept] >= std::numeric_limits<RealType>::epsilon() &&
              step_context.get_random_gen()() < prob[i_accept] * std::exp(log_gb[i_accept] - log_gf[i_accept]))
          {
            crowd.incAccept();
            twf_accept_list.push_back(crowd.get_walker_twfs()[i_accept]);
            elec_accept_list.push_back(crowd.get_walker_elecs()[i_accept]);
          }
          else
          {
            crowd.incReject();
            twf_reject_list.push_back(crowd.get_walker_twfs()[i_accept]);
            elec_reject_list.push_back(crowd.get_walker_elecs()[i_accept]);
          }

        TrialWaveFunction::flex_acceptMove(twf_accept_list, elec_accept_list, iat);
        TrialWaveFunction::flex_rejectMove(twf_reject_list, iat);

        ParticleSet::flex_acceptMove(elec_accept_list, iat);
        ParticleSet::flex_rejectMove(elec_reject_list, iat);
      }
    }
    std::for_each(crowd.get_walker_twfs().begin(), crowd.get_walker_twfs().end(),
                  [](TrialWaveFunction& twf) { twf.completeUpdates(); });
  }

  ParticleSet::flex_donePbyP(crowd.get_walker_elecs());
  timers.movepbyp_timer.stop();

  timers.buffer_timer.start();
  TrialWaveFunction::flex_updateBuffer(crowd.get_walker_twfs(), crowd.get_walker_elecs(), crowd.get_mcp_wfbuffers());

  auto saveElecPosAndGLToWalkers = [](ParticleSet& pset, ParticleSet::Walker_t& walker) { pset.saveWalker(walker); };
  for (int iw = 0; iw < crowd.size(); ++iw)
    saveElecPosAndGLToWalkers(walker_elecs[iw], walkers[iw]);
  timers.buffer_timer.stop();

  timers.hamiltonian_timer.start();
  auto& walker_hamiltonians = crowd.get_walker_hamiltonians();
  std::vector<QMCHamiltonian::FullPrecRealType> local_energies(
      QMCHamiltonian::flex_evaluate(walker_hamiltonians, walker_elecs));
  timers.hamiltonian_timer.stop();

  auto resetSigNLocalEnergy = [](MCPWalker& walker, TrialWaveFunction& twf, auto& local_energy) {
    walker.resetProperty(twf.getLogPsi(), twf.getPhase(), local_energy);
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    resetSigNLocalEnergy(walkers[iw], walker_twfs[iw], local_energies[iw]);

  // moved to be consistent with DMC
  timers.collectables_timer.start();
  auto evaluateNonPhysicalHamiltonianElements = [](QMCHamiltonian& ham, ParticleSet& pset, MCPWalker& walker) {
    ham.auxHevaluate(pset, walker);
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    evaluateNonPhysicalHamiltonianElements(walker_hamiltonians[iw], walker_elecs[iw], walkers[iw]);

  auto savePropertiesIntoWalker = [](QMCHamiltonian& ham, MCPWalker& walker) {
    ham.saveProperty(walker.getPropertyBase());
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    savePropertiesIntoWalker(walker_hamiltonians[iw], walkers[iw]);
  timers.collectables_timer.stop();
  // TODO:
  //  check if all moves failed
}


/** Thread body for VMC step
 *
 */
void VMCBatched::runVMCStep(int crowd_id,
                            const StateForThread& sft,
                            DriverTimers& timers,
                            std::vector<std::unique_ptr<ContextForSteps>>& context_for_steps,
                            std::vector<std::unique_ptr<Crowd>>& crowds)
{
  Crowd& crowd = *(crowds[crowd_id]);
  crowd.setRNGForHamiltonian(context_for_steps[crowd_id]->get_random_gen());

  int max_steps = sft.qmcdrv_input.get_max_steps();
  bool is_recompute_block =
      sft.recomputing_blocks ? (1 + sft.block) % sft.qmcdrv_input.get_blocks_between_recompute() == 0 : false;
  RealType cnorm = 1.0 / static_cast<RealType>(crowd.size());
  IndexType step = sft.step;
  // Are we entering the the last step of a block to recompute at?
  bool recompute_this_step = (is_recompute_block && (step + 1) == max_steps);
  advanceWalkers(sft, crowd, timers, *context_for_steps[crowd_id], recompute_this_step);
  crowd.accumulate(sft.population.get_num_global_walkers());
}

/** Runs the actual VMC section
 *
 *  Dependent on base class state machine
 *  Assumes state already updated from the following calls:
 *  1. QMCDriverNew::setStatus
 *  2. QMCDriverNew::putWalkers
 *  3. QMCDriverNew::process
 *
 *  At the moment I don't care about 1st touch, prove it matters
 *  If does consider giving more to the thread by value that should
 *  end up thread local. (I think)
 */
bool VMCBatched::run()
{
  IndexType num_blocks = qmcdriver_input_.get_max_blocks();
  //start the main estimator
  estimator_manager_->start(num_blocks);

  StateForThread vmc_state(qmcdriver_input_, vmcdriver_input_, *drift_modifier_, population_);

  LoopTimer vmc_loop;
  RunTimeControl runtimeControl(RunTimeManager, MaxCPUSecs);

  { // walker initialization
    ScopedTimer local_timer(&(timers_.init_walkers_timer));
    TasksOneToOne<> section_start_task(num_crowds_);
    section_start_task(initialLogEvaluation, std::ref(crowds_), std::ref(step_contexts_));
  }

  TasksOneToOne<> crowd_task(num_crowds_);

  auto runWarmupStep = [](int crowd_id, StateForThread& sft, DriverTimers& timers,
                          UPtrVector<ContextForSteps>& context_for_steps, UPtrVector<Crowd>& crowds) {
    Crowd& crowd = *(crowds[crowd_id]);
    advanceWalkers(sft, crowd, timers, *context_for_steps[crowd_id], false);
  };

  for (int step = 0; step < qmcdriver_input_.get_warmup_steps(); ++step)
  {
    ScopedTimer local_timer(&(timers_.run_steps_timer));
    crowd_task(runWarmupStep, vmc_state, std::ref(timers_), std::ref(step_contexts_), std::ref(crowds_));
  }

  for (int block = 0; block < num_blocks; ++block)
  {
    vmc_loop.start();
    vmc_state.recalculate_properties_period =
        (qmc_driver_mode_[QMC_UPDATE_MODE]) ? qmcdriver_input_.get_recalculate_properties_period() : 0;
    vmc_state.recomputing_blocks = qmcdriver_input_.get_blocks_between_recompute();

    estimator_manager_->startBlock(qmcdriver_input_.get_max_steps());

    for (auto& crowd : crowds_)
      crowd->startBlock(qmcdriver_input_.get_max_steps());
    for (int step = 0; step < qmcdriver_input_.get_max_steps(); ++step)
    {
      ScopedTimer local_timer(&(timers_.run_steps_timer));
      vmc_state.step = step;
      crowd_task(runVMCStep, vmc_state, timers_, std::ref(step_contexts_), std::ref(crowds_));
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

  // This is confusing logic from VMC.cpp want this functionality write documentation of this
  // and clean it up
  // bool wrotesamples = qmcdriver_input_.get_dump_config();
  // if (qmcdriver_input_.get_dump_config())
  // {
  //wrotesamples = W.dumpEnsemble(wClones, wOut, myComm->size(), nBlocks);
  //if (wrotesamples)
  //  app_log() << "  samples are written to the config.h5" << std::endl;
  // }

  // second argument was !wrotesample so if W.dumpEnsemble returns false or
  // dump_config is false from input then dump_walkers
  return finalize(num_blocks, true);
}


} // namespace qmcplusplus
