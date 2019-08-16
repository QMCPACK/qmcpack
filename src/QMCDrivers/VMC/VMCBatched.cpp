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

namespace qmcplusplus
{
/** Constructor maintains proper ownership of input parameters
   */
VMCBatched::VMCBatched(QMCDriverInput& qmcdriver_input,
                       VMCDriverInput& input,
                       MCPopulation& pop,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       WaveFunctionPool& ppool,
                       Communicate* comm)
    : QMCDriverNew(qmcdriver_input, pop, psi, h, ppool, comm), vmcdriver_input_(input)
{
  // qmc_driver_mode.set(QMC_UPDATE_MODE, 1);
  // qmc_driver_mode.set(QMC_WARMUP, 0);
}

VMCBatched::IndexType VMCBatched::calc_default_local_walkers()
{
  IndexType rw = vmcdriver_input_.get_requested_walkers_per_rank();
  if (rw < num_crowds_)
    rw = num_crowds_;
  walkers_per_crowd_      = (rw % num_crowds_) ? rw / num_crowds_ + 1 : rw / num_crowds_;
  IndexType local_walkers = walkers_per_crowd_ * num_crowds_;
  population_.set_num_local_walkers(local_walkers);
  population_.set_num_global_walkers(local_walkers * population_.get_num_ranks());
  if (rw != vmcdriver_input_.get_requested_walkers_per_rank())
      app_warning() << "VMCBatched driver has adjusted walkers per rank to: "
                    << local_walkers << '\n';

  if (vmcdriver_input_.get_samples() >= 0 ||
      vmcdriver_input_.get_samples_per_thread() >= 0 ||
      vmcdriver_input_.get_steps_between_samples() >= 0)
    app_warning() << "VMCBatched currently ignores samples and samplesperthread\n";

  // TODO: Simplify samples, samples per thread etc in the unified driver
  // see login in original VMC.cpp
  return local_walkers;
}

bool VMCBatched::run()
{
  //   //start the main estimator
  //   Estimators->start(nBlocks);
  //   for (int ip = 0; ip < NumThreads; ++ip)
  //     Movers[ip]->startRun(nBlocks, false);
  // #if !defined(REMOVE_TRACEMANAGER)
  //   Traces->startRun(nBlocks, traceClones);
  // #endif

  //   LoopTimer vmc_loop;
  //   RunTimeControl runtimeControl(RunTimeManager, MaxCPUSecs);
  //   bool enough_time_for_next_iteration = true;

  //   const bool has_collectables = W.Collectables.size();
  //   for (int block = 0; block < nBlocks; ++block)
  //   {
  //     vmc_loop.start();
  // #pragma omp parallel
  //     {
  //       int ip = omp_get_thread_num();
  //       //IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
  //       IndexType updatePeriod = (QMCDriverMode[QMC_UPDATE_MODE]) ? Period4CheckProperties : 0;

  //       //Is there any rebalance per block with VMC?  I don't think so.
  //       //Therefore there is no need to return to a global population
  //       MCWalkerConfiguration::iterator wit(W.begin() + wPerNode[ip]), wit_end(W.begin() + wPerNode[ip + 1]);
  //       Movers[ip]->startBlock(nSteps);
  //       int now_loc    = CurrentStep;
  //       RealType cnorm = 1.0 / static_cast<RealType>(wPerNode[ip + 1] - wPerNode[ip]);
  //       for (int step = 0; step < nSteps; ++step)
  //       {
  //         // Why?
  //         Movers[ip]->set_step(now_loc);
  //         //collectables are reset, it is accumulated while advancing walkers
  //         wClones[ip]->resetCollectables();
  //         bool recompute = (nBlocksBetweenRecompute && (step + 1) == nSteps &&
  //                           (1 + block) % nBlocksBetweenRecompute == 0 && QMCDriverMode[QMC_UPDATE_MODE]);
  //         Movers[ip]->advanceWalkers(wit, wit_end, recompute);
  //         if (has_collectables)
  //           wClones[ip]->Collectables *= cnorm;
  //         Movers[ip]->accumulate(wit, wit_end);
  //         ++now_loc;
  //         if (Period4WalkerDump && now_loc % Period4WalkerDump == 0)
  //           wClones[ip]->saveEnsemble(wit, wit_end);
  //         //           if(storeConfigs && (now_loc%storeConfigs == 0))
  //         //             ForwardWalkingHistory.storeConfigsForForwardWalking(*wClones[ip]);
  //       }
  //       Movers[ip]->stopBlock(false);
  //     } //end-of-parallel for
  //     //Estimators->accumulateCollectables(wClones,nSteps);
  //     CurrentStep += nSteps;
  //     Estimators->stopBlock(estimatorClones);
  // #if !defined(REMOVE_TRACEMANAGER)
  //     Traces->write_buffers(traceClones, block);
  // #endif
  //     if (storeConfigs)
  //       recordBlock(block);
  //     vmc_loop.stop();
  //     enough_time_for_next_iteration = runtimeControl.enough_time_for_next_iteration(vmc_loop);
  //     myComm->bcast(enough_time_for_next_iteration);
  //     if (!enough_time_for_next_iteration)
  //     {
  //       app_log() << runtimeControl.time_limit_message("VMC", block);
  //       break;
  //     }
  //   } //block
  //   Estimators->stop(estimatorClones);
  //   for (int ip = 0; ip < NumThreads; ++ip)
  //     Movers[ip]->stopRun2();
  // #if !defined(REMOVE_TRACEMANAGER)
  //   Traces->stopRun();
  // #endif
  //   //copy back the random states
  // #ifndef USE_FAKE_RNG
  //   for (int ip = 0; ip < NumThreads; ++ip)
  //     *(RandomNumberControl::Children[ip]) = *(Rng[ip]);
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
