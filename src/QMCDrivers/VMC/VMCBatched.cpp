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

void VMCBatched::runVMCBlock(int crowd_id, const QMCDriverInput& qmcdriver_input, std::vector<Crowd>& crowds)
{  
  int nAccept              = 0;
  int nReject              = 0;
  int nAllRejected         = 0;
  int nNodeCrossing        = 0;
  int NonLocalMoveAccepted = 0;

  crowds[crowd_id].startBlock(qmcdriver_input.get_max_steps());
  //Is there any rebalance per block with VMC?  I don't think so.
  //Therefore there is no need to return to a global population
  // MCWalkerConfiguration::iterator wit(W.begin() + wPerNode[crowd_id]), wit_end(W.begin() + wPerNode[crowd_id + 1]);
//     Movers[crowd_id]->startBlock(nSteps);
//     int now_loc    = CurrentStep;
//     RealType cnorm = 1.0 / static_cast<RealType>(wPerNode[crowd_id + 1] - wPerNode[crowd_id]);
//     for (int step = 0; step < nSteps; ++step)
//     {
//       // Why?
//       Movers[crowd_id]->set_step(now_loc);
//       //collectables are reset, it is accumulated while advancing walkers
//       wClones[crowd_id]->resetCollectables();
//       bool recompute = (nBlocksBetweenRecompute && (step + 1) == nSteps && (1 + block) % nBlocksBetweenRecompute == 0 &&
//                         QMCDriverMode[QMC_UPDATE_MODE]);
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
 */
bool VMCBatched::run()
{
  IndexType num_blocks = qmcdriver_input_.get_max_blocks();
  //start the main estimator
  estimator_manager_->start(num_blocks);

  LoopTimer vmc_loop;
  RunTimeControl runtimeControl(RunTimeManager, MaxCPUSecs);

  // TODO: Do collectables need to be rethought
  //   const bool has_collectables = W.Collectables.size();

  for (int block = 0; block < num_blocks; ++block)
  {
    vmc_loop.start();
    IndexType updatePeriod = (qmc_driver_mode_[QMC_UPDATE_MODE]) ? qmcdriver_input_.get_recalculate_properties_period() : 0;

    estimator_manager_->startBlock(qmcdriver_input_.get_max_steps());

    TasksOneToOne<> crowd_task(num_crowds_);
    crowd_task(runVMCBlock, qmcdriver_input_, crowds_);
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
