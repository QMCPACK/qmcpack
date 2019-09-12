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
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "Concurrency/TasksOneToOne.hpp"
#include "Concurrency/Info.hpp"
#include "Utilities/RunTimeManager.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include <boost/range/combine.hpp>

namespace qmcplusplus
{
/** Constructor maintains proper ownership of input parameters
   */
VMCBatched::VMCBatched(QMCDriverInput&& qmcdriver_input,
                       VMCDriverInput&& input,
                       MCPopulation&& pop,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       WaveFunctionPool& ppool,
                       Communicate* comm)
    : QMCDriverNew(std::move(qmcdriver_input), std::move(pop), psi, h, ppool, comm), vmcdriver_input_(input)
{
  // qmc_driver_mode.set(QMC_UPDATE_MODE, 1);
  // qmc_driver_mode.set(QMC_WARMUP, 0);
}

VMCBatched::IndexType VMCBatched::calc_default_local_walkers()
{
  int num_threads(Concurrency::maxThreads<>());

  // Do to a work-around currently in QMCDriverNew::QMCDriverNew this should never be true.
  // I'm leaving this because this is what should happen for vmc.
  if (num_crowds_ > num_threads)
  {
    std::stringstream error_msg;
    error_msg << "Bad Input: num_crowds (" << qmcdriver_input_.get_num_crowds() << ") > num_threads (" << num_threads
              << ")\n";
    throw std::runtime_error(error_msg.str());
  }

  IndexType rw = vmcdriver_input_.get_requested_walkers_per_rank();
  if (rw < num_crowds_)
    rw = num_crowds_;
  walkers_per_crowd_      = (rw % num_crowds_) ? rw / num_crowds_ + 1 : rw / num_crowds_;
  IndexType local_walkers = walkers_per_crowd_ * num_crowds_;
  population_.set_num_local_walkers(local_walkers);
  population_.set_num_global_walkers(local_walkers * population_.get_num_ranks());
  if (rw != vmcdriver_input_.get_requested_walkers_per_rank())
    app_warning() << "VMCBatched driver has adjusted walkers per rank to: " << local_walkers << '\n';

  if (vmcdriver_input_.get_samples() >= 0 || vmcdriver_input_.get_samples_per_thread() >= 0 ||
      vmcdriver_input_.get_steps_between_samples() >= 0)
    app_warning() << "VMCBatched currently ignores samples and samplesperthread\n";

  // TODO: Simplify samples, samples per thread etc in the unified driver
  // see login in original VMC.cpp
  return local_walkers;
}

void VMCBatched::advanceWalkers(const StateForThread& sft, Crowd& crowd, ContextForSteps& step_context, bool recompute)
{
  crowd.loadWalkers();
  auto it_walker_twfs  = crowd.beginTrialWaveFunctions();
  auto it_mcp_walkers  = crowd.beginWalkers();
  auto it_walker_elecs = crowd.beginElectrons();
  while (it_walker_twfs != crowd.endTrialWaveFunctions())
  {
    it_walker_twfs->get().copyFromBuffer(it_walker_elecs->get(), it_mcp_walkers->get().DataSet);
  }

  int num_walkers = crowd.size();
  int num_particles = sft.population.get_num_particles();
  // Note std::vector<bool> is not like the rest of stl.
  std::vector<bool> moved(num_walkers, false);
  constexpr RealType mhalf(-0.5);
  bool use_drift = sft.vmcdrv_input.get_use_drift();
  //This generates an entire steps worth of deltas.
  step_context.nextDeltaRs();
  auto it_delta_r = step_context.deltaRsBegin();
  std::vector<PosType> drifts(num_walkers);

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
      crowd.clearResults();
      ParticleSet::flex_setActive(crowd.get_walker_elecs(), iat);
      // step_context.deltaRsBegin returns an iterator to a flat series of PosTypes
      // fastest in walkers then particles
      auto delta_r_start = it_delta_r + iat * num_particles * num_walkers;

      if (use_drift)
      {
        TrialWaveFunction::flex_evalGrad(crowd.get_walker_twfs(), crowd.get_walker_elecs(), iat, crowd.get_grads_now());
        sft.drift_modifier.getDrifts(tauovermass, crowd.get_grads_now(), drifts);
        std::transform(drifts.begin(),drifts.end(),
                       delta_r_start, drifts.begin(),
                       [sqrttau](PosType& drift, PosType& delta_r){
                         return drift + (sqrttau * delta_r);});
      }
      else
      {
        std::transform(drifts.begin(),drifts.end(),delta_r_start, drifts.begin(),
                       [sqrttau](auto& drift, auto& delta_r){
                         return sqrttau * delta_r;});
      }

      auto elecs = crowd.get_walker_elecs();
      ParticleSet::flex_makeMove(crowd.get_walker_elecs(), iat, drifts);

      // This is inelegant
      if(use_drift)
      {
        TrialWaveFunction::flex_ratioGrad(crowd.get_walker_twfs(),crowd.get_walker_elecs(), iat, crowd.get_ratios(), crowd.get_grads_new());
        auto delta_r_end = it_delta_r + iat * num_particles * num_walkers + num_walkers;
        std::transform(delta_r_start,
                       delta_r_end,
                       crowd.get_log_gf().begin(),
                       [mhalf](auto& delta_r){ return mhalf * dot(delta_r, delta_r); });
        sft.drift_modifier.getDrifts(sft.qmcdrv_input.get_tau(), crowd.get_grads_new(), drifts);
        
        std::transform(crowd.beginElectrons(), crowd.endElectrons(),
                       drifts.begin(), drifts.begin(),
                       [iat](auto& elecs, auto& drift){
                         return elecs.get().R[iat] - elecs.get().activePos - drift; });

        std::transform(drifts.begin(),
                       drifts.end(),
                       crowd.get_log_gb().begin(),
                       [oneover2tau](auto& drift){ return - oneover2tau * dot(drift, drift);});        
      }
      else
      {
        TrialWaveFunction::flex_calcRatio(crowd.get_walker_twfs(), crowd.get_walker_elecs(),
                                      iat, crowd.get_ratios());
      }

      std::transform(crowd.get_ratios().begin(),
                     crowd.get_ratios().end(),
                     crowd.get_prob().begin(),
                     [](auto ratio){ return std::real(ratio) * std::real(ratio); });

      for(int i_accept = 0; i_accept < num_walkers; ++i_accept)
      {
        TrialWaveFunction& walker_twf = crowd.get_walker_twfs()[i_accept].get();
        ParticleSet& walker_elecs = crowd.get_walker_elecs()[i_accept].get();
        auto prob = crowd.get_prob()[i_accept];
        auto log_gf = crowd.get_log_gf()[i_accept];
        auto log_gb = crowd.get_log_gb()[i_accept];
          
        if (prob >= std::numeric_limits<RealType>::epsilon()
            && step_context.get_random_gen()() < prob * std::exp(log_gf - log_gb))
        {
          step_context.incAccept();
          walker_twf.acceptMove(walker_elecs, iat);
          walker_elecs.acceptMove(iat);
        }
        else
        {
          step_context.incReject();
          walker_elecs.rejectMove(iat);
          walker_twf.rejectMove(iat);
        }
      }
    }
  }
  std::for_each(crowd.get_walker_twfs().begin(),
                crowd.get_walker_twfs().end(),
                [](auto& twf){twf.get().completeUpdates();});
  std::for_each(crowd.get_walker_elecs().begin(),
                crowd.get_walker_elecs().end(),
                [](auto& elecs){elecs.get().donePbyP();});
  // TODO:
  //  update the buffer.
  //  save the walkers
  //  evaluate the Hamiltonian
  //  accumulate
  //  check if all moves failed

}

/** Thread body for VMC block
 *
 *  Things to consider:
 *  - should qmcdriver_input be a copy local to the core in Crowd
 */
    void VMCBatched::runVMCStep(int crowd_id, const StateForThread& sft,
                                std::vector<std::unique_ptr<ContextForSteps>>& context_for_steps,
                                std::vector<std::unique_ptr<Crowd>>& crowds)
    {
      int nAccept              = 0;
      int nReject              = 0;
      int nAllRejected         = 0;
      int nNodeCrossing        = 0;
      int NonLocalMoveAccepted = 0;

      Crowd& crowd = *(crowds[crowd_id]);
      crowd.startBlock(sft.qmcdrv_input.get_max_steps());

      int max_steps = sft.qmcdrv_input.get_max_steps();
      bool is_recompute_block =
          sft.recomputing_blocks ? (1 + sft.block) % sft.qmcdrv_input.get_blocks_between_recompute() == 0 : false;
      RealType cnorm = 1.0 / static_cast<RealType>(crowd.size());
      IndexType step = sft.step;
      // Are we entering the the last step of a block to recompute at?
      bool recompute_this_step = (is_recompute_block && (step + 1) == max_steps);
      advanceWalkers(sft, crowd, *context_for_steps[crowd_id], recompute_this_step);
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

      // TODO: Do collectables need to be rethought
      //   const bool has_collectables = W.Collectables.size();

      for (int block = 0; block < num_blocks; ++block)
      {
        vmc_loop.start();
        vmc_state.recalculate_properties_period =
            (qmc_driver_mode_[QMC_UPDATE_MODE]) ? qmcdriver_input_.get_recalculate_properties_period() : 0;
        vmc_state.recomputing_blocks = qmcdriver_input_.get_blocks_between_recompute();

        estimator_manager_->startBlock(qmcdriver_input_.get_max_steps());

        for (int step = 0; step < qmcdriver_input_.get_max_steps(); ++step)
        {
          vmc_state.step = step;
          TasksOneToOne<> crowd_task(num_crowds_);
          crowd_task(runVMCStep, vmc_state, std::ref(step_contexts_), std::ref(crowds_));
        }
      }

      return false;
    }


  } // namespace qmcplusplus
