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
    : QMCDriverNew(std::move(qmcdriver_input), pop, psi, h, ppool, comm), vmcdriver_input_(input)
{
  // qmc_driver_mode.set(QMC_UPDATE_MODE, 1);
  // qmc_driver_mode.set(QMC_WARMUP, 0);
}

VMCBatched::IndexType VMCBatched::calc_default_local_walkers()
{
  int num_threads(Concurrency::maxThreads<>());

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

void VMCBatched::advanceWalkers(Crowd& crowd, MoveContext& move_context, bool recompute)
{
  move_context.loadCrowd(crowd);
//   Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
//   Psi.copyFromBuffer(W, w_buffer);
//   myTimers[0]->stop();
  
//   // start PbyP moves
//   myTimers[1]->start();

  bool moved = false;
  constexpr RealType mhalf(-0.5);
//   for (int iter = 0; iter < nSubSteps; ++iter)
//   {
//     //create a 3N-Dimensional Gaussian with variance=1


  // up and down electrons are "species" within qmpack
  for (int ig = 0; ig < move_context.get_num_groups(); ++ig) //loop over species
  {}
//       RealType tauovermass = Tau * MassInvS[ig];
//       RealType oneover2tau = 0.5 / (tauovermass);
//       RealType sqrttau     = std::sqrt(tauovermass);
//       for (int iat = W.first(ig); iat < W.last(ig); ++iat)
//       {
//         W.setActive(iat);
//         PosType dr;
//         if (UseDrift)
//         {
//           GradType grad_now = Psi.evalGrad(W, iat);
//           DriftModifier->getDrift(tauovermass, grad_now, dr);
//           dr += sqrttau * deltaR[iat];
//         }
//         else
//         {
//           dr = sqrttau * deltaR[iat];
//         }
//         if (!W.makeMoveAndCheck(iat, dr))
//         {
//           ++nReject;
//           continue;
//         }
//         RealType logGf(1), logGb(1), prob;
//         if (UseDrift)
//         {
//           GradType grad_new;
//           RealType ratio = Psi.ratioGrad(W, iat, grad_new);
//           prob           = ratio * ratio;
//           logGf          = mhalf * dot(deltaR[iat], deltaR[iat]);
//           DriftModifier->getDrift(tauovermass, grad_new, dr);
//           dr    = W.R[iat] - W.activePos - dr;
//           logGb = -oneover2tau * dot(dr, dr);
//         }
//         else
//         {
//           RealType ratio = Psi.ratio(W, iat);
//           prob           = ratio * ratio;
//         }
//         if (prob >= std::numeric_limits<RealType>::epsilon() && RandomGen() < prob * std::exp(logGb - logGf))
//         {
//           moved = true;
//           ++nAccept;
//           Psi.acceptMove(W, iat);
//           W.acceptMove(iat);
//         }
//         else
//         {
//           ++nReject;
//           W.rejectMove(iat);
//           Psi.rejectMove(iat);
//         }
//       }
//     }
//     Psi.completeUpdates();
//   }
//   W.donePbyP();
//   myTimers[1]->stop();
//   myTimers[0]->start();
//   RealType logpsi = Psi.updateBuffer(W, w_buffer, recompute);
//   W.saveWalker(thisWalker);
//   myTimers[0]->stop();
//   // end PbyP moves
//   myTimers[2]->start();
//   FullPrecRealType eloc = H.evaluate(W);
//   thisWalker.resetProperty(logpsi, Psi.getPhase(), eloc);
//   myTimers[2]->stop();
//   myTimers[3]->start();
//   H.auxHevaluate(W, thisWalker);
//   H.saveProperty(thisWalker.getPropertyBase());
//   myTimers[3]->stop();
// #if !defined(REMOVE_TRACEMANAGER)
//   Traces->buffer_sample(W.current_step);
// #endif
//   if (!moved)
//     ++nAllRejected;
  
}

/** Thread body for VMC block
 *
 *  Things to consider:
 *  - should qmcdriver_input be a copy local to the core in Crowd
 */
void VMCBatched::runVMCStep(int crowd_id,
                            const QMCDriverInput& qmcdriver_input,
                            const StateForThreads vmc_state,
                            std::vector<std::unique_ptr<MoveContext>>& move_context,
                            std::vector<Crowd>& crowds)
{
  int nAccept              = 0;
  int nReject              = 0;
  int nAllRejected         = 0;
  int nNodeCrossing        = 0;
  int NonLocalMoveAccepted = 0;

  Crowd& crowd = crowds[crowd_id];
  crowd.startBlock(qmcdriver_input.get_max_steps());
  
  int max_steps = qmcdriver_input.get_max_steps();
  bool is_recompute_block = vmc_state.recomputing_blocks ? (1 + vmc_state.block) % qmcdriver_input.get_blocks_between_recompute() == 0 : false;
  RealType cnorm = 1.0 / static_cast<RealType>(crowd.size());
  IndexType step = vmc_state.step;
  // Are we entering the the last step of a block to recompute at?
  bool recompute_this_step = (is_recompute_block && (step + 1) == max_steps );
  advanceWalkers(crowd, *move_context[crowd_id], recompute_this_step);
  
  //       Movers[crowd_id]->advanceWalkers(wit, wit_end, recompute);
  //       if (has_collectables)
  //         wClones[crowd_id]->Collectables *= cnorm;
  //       Movers[crowd_id]->accumulate(wit, wit_end);
  //       ++now_loc;
  //       if (Period4WalkerDump && now_loc % Period4WalkerDump == 0)
  //         wClones[crowd_id]->saveEnsemble(wit, wit_end);
  //       //           if(storeConfigs && (now_loc%storeConfigs == 0))
  //       //             ForwardWalkingHistory.storeConfigsForForwardWalking(*wClones[crowd_id]);
  //     }
  //     Movers[crowd_id]->stopBlock(false);
  //   } //end-of-parallel for
  //   //Estimators->accumulateCollectables(wClones,nSteps);
  //   CurrentStep += nSteps;
  //   Estimators->stopBlock(estimatorClones);
  // #if !defined(REMOVE_TRACEMANAGER)
  //   Traces->write_buffers(traceClones, block);
  // #endif
  //   if (storeConfigs)
  //     recordBlock(block);
  //   vmc_loop.stop();
  //   enough_time_for_next_iteration = runtimeControl.enough_time_for_next_iteration(vmc_loop);
  //   myComm->bcast(enough_time_for_next_iteration);
  //   if (!enough_time_for_next_iteration)
  //   {
  //     app_log() << runtimeControl.time_limit_message("VMC", block);
  //     break;
  //   }
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

  StateForThreads vmc_state;

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

    for(int step = 0; step < qmcdriver_input_.get_max_steps(); ++step)
    {
      vmc_state.step =  step;
      TasksOneToOne<> crowd_task(num_crowds_);
      crowd_task(runVMCStep, std::cref(qmcdriver_input_), vmc_state, std::ref(move_contexts_), std::ref(crowds_));
    }
  }

  //   } //block
  //   Estimators->stop(estimatorClones);
  //   for (int crowd_id = 0; crowd_id < NumThreads; ++crowd_id)
  //     Movers[crowd_id]->stopRun2();
  // #if !defined(REMOVE_TRACEMANAGER)
  //   Traces->stopRun();
  // #endif
  //   //copy back the random states
  // #ifndef USE_FAKE_RNG
  //   for (int crowd_id = 0; crowd_id < NumThreads; ++crowd_id)
  //     *(RandomNumberControl::Children[crowd_id]) = *(Rng[crowd_id]);
  // #endif
  //   ///write samples to a file
  //   bool wrotesamples = DumpConfig;
  //   if (DumpConfig)
  //   {
  //     wrotesamples = W.dumpEnsemble(wClones, wOut, myComm->size(), nBlocks);
  //     if (wrotesamples)
  //       app_log() << "  samples are written to the config.h5" << std::endl;
  //   }
  //   //finalize a qmc section
  //   return finalize(nBlocks, !wrotesamples);
  // }
  return false;
}


} // namespace qmcplusplus
