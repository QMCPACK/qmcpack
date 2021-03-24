//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: VMC.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "VMCBatched.h"
#include "Concurrency/ParallelExecutor.hpp"
#include "Concurrency/Info.hpp"
#include "Message/UniformCommunicateError.h"
#include "Message/CommOperators.h"
#include "Utilities/RunTimeManager.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Particle/MCSample.h"
#include "MemoryUsage.h"

namespace qmcplusplus
{
/** Constructor maintains proper ownership of input parameters
   */
VMCBatched::VMCBatched(const ProjectData& project_data,
                       QMCDriverInput&& qmcdriver_input,
                       VMCDriverInput&& input,
                       MCPopulation&& pop,
                       SampleStack& samples,
                       Communicate* comm)
    : QMCDriverNew(project_data, std::move(qmcdriver_input), std::move(pop), "VMCBatched::", comm, "VMCBatched"),
      vmcdriver_input_(input),
      samples_(samples),
      collect_samples_(false)
{}

void VMCBatched::advanceWalkers(const StateForThread& sft,
                                Crowd& crowd,
                                QMCDriverNew::DriverTimers& timers,
                                ContextForSteps& step_context,
                                bool recompute)
{
  if (crowd.size() == 0)
    return;
  auto& ps_dispatcher  = crowd.dispatchers_.ps_dispatcher_;
  auto& twf_dispatcher = crowd.dispatchers_.twf_dispatcher_;
  auto& ham_dispatcher = crowd.dispatchers_.ham_dispatcher_;
  auto& walkers        = crowd.get_walkers();
  const RefVectorWithLeader<ParticleSet> walker_elecs(crowd.get_walker_elecs()[0], crowd.get_walker_elecs());
  const RefVectorWithLeader<TrialWaveFunction> walker_twfs(crowd.get_walker_twfs()[0], crowd.get_walker_twfs());

  ResourceCollectionTeamLock<ParticleSet> pset_res_lock(crowd.getSharedResource().pset_res, walker_elecs);
  DriverWalkerResourceCollectionLock pbyp_lock(crowd.getSharedResource(), crowd.get_walker_twfs()[0],
                                               crowd.get_walker_hamiltonians()[0]);

  assert(QMCDriverNew::checkLogAndGL(crowd));

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
  std::vector<bool> isAccepted;
  std::vector<std::reference_wrapper<TrialWaveFunction>> twf_accept_list, twf_reject_list;
  isAccepted.reserve(num_walkers);

  for (int sub_step = 0; sub_step < sft.qmcdrv_input.get_sub_steps(); sub_step++)
  {
    //This generates an entire steps worth of deltas.
    step_context.nextDeltaRs(num_walkers * sft.population.get_num_particles());

    // up and down electrons are "species" within qmpack
    for (int ig = 0; ig < step_context.get_num_groups(); ++ig) //loop over species
    {
      RealType tauovermass = sft.qmcdrv_input.get_tau() * sft.population.get_ptclgrp_inv_mass()[ig];
      RealType oneover2tau = 0.5 / (tauovermass);
      RealType sqrttau     = std::sqrt(tauovermass);

      twf_dispatcher.flex_prepareGroup(walker_twfs, walker_elecs, ig);

      int start_index = step_context.getPtclGroupStart(ig);
      int end_index   = step_context.getPtclGroupEnd(ig);
      for (int iat = start_index; iat < end_index; ++iat)
      {
        // step_context.deltaRsBegin returns an iterator to a flat series of PosTypes
        // fastest in walkers then particles
        auto delta_r_start = step_context.deltaRsBegin() + iat * num_walkers;
        auto delta_r_end   = delta_r_start + num_walkers;

        if (use_drift)
        {
          twf_dispatcher.flex_evalGrad(walker_twfs, walker_elecs, iat, grads_now);
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

        ps_dispatcher.flex_makeMove(walker_elecs, iat, drifts);

        // This is inelegant
        if (use_drift)
        {
          twf_dispatcher.flex_calcRatioGrad(walker_twfs, walker_elecs, iat, ratios, grads_new);
          std::transform(delta_r_start, delta_r_end, log_gf.begin(),
                         [](const PosType& delta_r) { return mhalf * dot(delta_r, delta_r); });

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
          twf_dispatcher.flex_calcRatio(walker_twfs, walker_elecs, iat, ratios);
        }

        std::transform(ratios.begin(), ratios.end(), prob.begin(), [](auto ratio) { return std::norm(ratio); });

        isAccepted.clear();

        for (int i_accept = 0; i_accept < num_walkers; ++i_accept)
          if (prob[i_accept] >= std::numeric_limits<RealType>::epsilon() &&
              step_context.get_random_gen()() < prob[i_accept] * std::exp(log_gb[i_accept] - log_gf[i_accept]))
          {
            crowd.incAccept();
            isAccepted.push_back(true);
          }
          else
          {
            crowd.incReject();
            isAccepted.push_back(false);
          }

        twf_dispatcher.flex_accept_rejectMove(walker_twfs, walker_elecs, iat, isAccepted, true);

        ps_dispatcher.flex_accept_rejectMove(walker_elecs, iat, isAccepted);
      }
    }
    twf_dispatcher.flex_completeUpdates(walker_twfs);
  }

  ps_dispatcher.flex_donePbyP(walker_elecs);
  timers.movepbyp_timer.stop();

  timers.buffer_timer.start();
  twf_dispatcher.flex_evaluateGL(walker_twfs, walker_elecs, recompute);

  assert(QMCDriverNew::checkLogAndGL(crowd));
  timers.buffer_timer.stop();

  timers.hamiltonian_timer.start();
  const RefVectorWithLeader<QMCHamiltonian> walker_hamiltonians(crowd.get_walker_hamiltonians()[0],
                                                                crowd.get_walker_hamiltonians());
  std::vector<QMCHamiltonian::FullPrecRealType> local_energies(
      ham_dispatcher.flex_evaluate(walker_hamiltonians, walker_elecs));
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
  // \todo delete
  RealType cnorm = 1.0 / static_cast<RealType>(crowd.size());
  IndexType step = sft.step;
  // Are we entering the the last step of a block to recompute at?
  bool recompute_this_step = (sft.is_recomputing_block && (step + 1) == max_steps);
  advanceWalkers(sft, crowd, timers, *context_for_steps[crowd_id], recompute_this_step);
  crowd.accumulate();
}

void VMCBatched::process(xmlNodePtr node)
{
  print_mem("VMCBatched before initialization", app_log());
  // \todo get total walkers should be coming from VMCDriverInpu
  try
  {
    QMCDriverNew::AdjustedWalkerCounts awc =
        adjustGlobalWalkerCount(myComm->size(), myComm->rank(), qmcdriver_input_.get_total_walkers(),
                                qmcdriver_input_.get_walkers_per_rank(), 1.0, qmcdriver_input_.get_num_crowds());

    Base::startup(node, awc);
  }
  catch (const UniformCommunicateError& ue)
  {
    myComm->barrier_and_abort(ue.what());
  }
}

int VMCBatched::compute_samples_per_rank(const QMCDriverInput& qmcdriver_input, const IndexType local_walkers)
{
  int nblocks = qmcdriver_input.get_max_blocks();
  int nsteps  = qmcdriver_input.get_max_steps();
  return nblocks * nsteps * local_walkers;
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
  estimator_manager_->startDriverRun();

  StateForThread vmc_state(qmcdriver_input_, vmcdriver_input_, *drift_modifier_, population_);

  LoopTimer<> vmc_loop;
  RunTimeControl<> runtimeControl(run_time_manager, project_data_.getMaxCPUSeconds(), project_data_.getTitle(),
                                  myComm->rank() == 0);

  { // walker initialization
    ScopedTimer local_timer(timers_.init_walkers_timer);
    ParallelExecutor<> section_start_task;
    section_start_task(crowds_.size(), initialLogEvaluation, std::ref(crowds_), std::ref(step_contexts_));
  }

  print_mem("VMCBatched after initialLogEvaluation", app_summary());

  ParallelExecutor<> crowd_task;

  if (qmcdriver_input_.get_warmup_steps() > 0)
  {
    // Run warm-up steps
    auto runWarmupStep = [](int crowd_id, StateForThread& sft, DriverTimers& timers,
                            UPtrVector<ContextForSteps>& context_for_steps, UPtrVector<Crowd>& crowds) {
      Crowd& crowd = *(crowds[crowd_id]);
      advanceWalkers(sft, crowd, timers, *context_for_steps[crowd_id], false);
    };

    for (int step = 0; step < qmcdriver_input_.get_warmup_steps(); ++step)
    {
      ScopedTimer local_timer(timers_.run_steps_timer);
      crowd_task(crowds_.size(), runWarmupStep, vmc_state, std::ref(timers_), std::ref(step_contexts_),
                 std::ref(crowds_));
    }

    app_log() << "Warm-up is completed!" << std::endl;
    print_mem("VMCBatched after Warmup", app_log());
  }

  for (int block = 0; block < num_blocks; ++block)
  {
    vmc_loop.start();
    vmc_state.recalculate_properties_period =
        (qmc_driver_mode_[QMC_UPDATE_MODE]) ? qmcdriver_input_.get_recalculate_properties_period() : 0;
    vmc_state.is_recomputing_block = qmcdriver_input_.get_blocks_between_recompute()
        ? (1 + block) % qmcdriver_input_.get_blocks_between_recompute() == 0
        : false;

    estimator_manager_->startBlock(qmcdriver_input_.get_max_steps());

    for (auto& crowd : crowds_)
      crowd->startBlock(qmcdriver_input_.get_max_steps());
    for (int step = 0; step < qmcdriver_input_.get_max_steps(); ++step)
    {
      ScopedTimer local_timer(timers_.run_steps_timer);
      vmc_state.step = step;
      crowd_task(crowds_.size(), runVMCStep, vmc_state, timers_, std::ref(step_contexts_), std::ref(crowds_));

      if (collect_samples_)
      {
        const auto& elec_psets = population_.get_elec_particle_sets();
        for (const auto& walker : elec_psets)
        {
          samples_.appendSample(MCSample(*walker));
        }
      }
    }
    print_mem("VMCBatched after a block", app_debug_stream());
    endBlock();
    vmc_loop.stop();

    bool stop_requested = false;
    // Rank 0 decides whether the time limit was reached
    if (!myComm->rank())
      stop_requested = runtimeControl.checkStop(vmc_loop);
    myComm->bcast(stop_requested);

    if (stop_requested)
    {
      if (!myComm->rank())
        app_log() << runtimeControl.generateStopMessage("VMCBatched", block);
      run_time_manager.markStop();
      break;
    }
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
  {
    std::ostringstream o;
    FullPrecRealType ene, var;
    estimator_manager_->getApproximateEnergyVariance(ene, var);
    o << "====================================================";
    o << "\n  End of a VMC block";
    o << "\n    QMC counter        = " << project_data_.getSeriesIndex();
    o << "\n    time step          = " << qmcdriver_input_.get_tau();
    o << "\n    reference energy   = " << ene;
    o << "\n    reference variance = " << var;
    o << "\n====================================================";
    app_log() << o.str() << std::endl;
  }

  print_mem("VMCBatched ends", app_log());

  estimator_manager_->stopDriverRun();

  return finalize(num_blocks, true);
}

void VMCBatched::enable_sample_collection()
{
  int samples = compute_samples_per_rank(qmcdriver_input_, population_.get_num_local_walkers());
  samples_.setMaxSamples(samples);
  collect_samples_ = true;

  int total_samples = samples * population_.get_num_ranks();
  app_log() << "VMCBatched Driver collecting samples, samples per rank = " << samples << '\n';
  app_log() << "                                      total samples    = " << total_samples << '\n';
}

} // namespace qmcplusplus
