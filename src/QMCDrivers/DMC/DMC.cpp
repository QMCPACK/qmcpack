//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/DMC/DMC.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Message/OpenMP.h"
#include "Utilities/Timer.h"
#include "Utilities/RunTimeManager.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Utilities/ProgressReportEngine.h"
#include <qmc_common.h>
#include "ADIOS/ADIOS_profile.h"
#include "Utilities/FairDivide.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif

namespace qmcplusplus
{
/// Constructor.
DMC::DMC(MCWalkerConfiguration& w,
         TrialWaveFunction& psi,
         QMCHamiltonian& h,
         WaveFunctionPool& ppool,
         Communicate* comm)
    : QMCDriver(w, psi, h, ppool, comm),
      KillNodeCrossing(0),
      Reconfiguration("no"),
      BranchInterval(-1),
      mover_MaxAge(-1)
{
  RootName = "dmc";
  QMCType  = "DMC";
  qmc_driver_mode.set(QMC_UPDATE_MODE, 1);
  m_param.add(KillWalker, "killnode", "string");
  m_param.add(Reconfiguration, "reconfiguration", "string");
  //m_param.add(BranchInterval,"branchInterval","string");
  m_param.add(NonLocalMove, "nonlocalmove", "string");
  m_param.add(NonLocalMove, "nonlocalmoves", "string");
  m_param.add(mover_MaxAge, "MaxAge", "double");
  //DMC overwrites ConstPopulation
  ConstPopulation = false;
}

void DMC::resetUpdateEngines()
{
  ReportEngine PRE("DMC", "resetUpdateEngines");
  bool fixW = (Reconfiguration == "yes");
  makeClones(W, Psi, H);
  Timer init_timer;
  if (Movers.empty())
  {
    W.loadEnsemble(wClones);
    setWalkerOffsets();
    int nw_multi = branchEngine->initWalkerController(W, fixW, false);
    if (nw_multi > 1)
    {
      W.createWalkers((nw_multi - 1) * W.getActiveWalkers());
      setWalkerOffsets();
    }
    //if(qmc_driver_mode[QMC_UPDATE_MODE]) W.clearAuxDataSet();
    Movers.resize(NumThreads, 0);
    Rng.resize(NumThreads, 0);
    estimatorClones.resize(NumThreads, 0);
    traceClones.resize(NumThreads, 0);
    FairDivideLow(W.getActiveWalkers(), NumThreads, wPerNode);
    {
      //log file
      std::ostringstream o;
      o << "  Initial partition of walkers on a node: ";
      copy(wPerNode.begin(), wPerNode.end(), std::ostream_iterator<int>(o, " "));
      o << "\n";
      if (qmc_driver_mode[QMC_UPDATE_MODE])
        o << "  Updates by particle-by-particle moves";
      else
        o << "  Updates by walker moves";
      if (KillNodeCrossing)
        o << "\n  Walkers are killed when a node crossing is detected";
      else
        o << "\n  DMC moves are rejected when a node crossing is detected";
      app_log() << o.str() << std::endl;
    }
#pragma omp parallel for
    for (int ip = 0; ip < NumThreads; ++ip)
    {
      estimatorClones[ip] = new EstimatorManagerBase(*Estimators);
      estimatorClones[ip]->setCollectionMode(false);
#if !defined(REMOVE_TRACEMANAGER)
      traceClones[ip] = Traces->makeClone();
#endif
#ifdef USE_FAKE_RNG
      Rng[ip] = new FakeRandom();
#else
      Rng[ip] = new RandomGenerator_t(*RandomNumberControl::Children[ip]);
      hClones[ip]->setRandomGenerator(Rng[ip]);
#endif
      if (qmc_driver_mode[QMC_UPDATE_MODE])
      {
        Movers[ip] = new DMCUpdatePbyPWithRejectionFast(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip]);
        Movers[ip]->put(qmcNode);
        Movers[ip]->resetRun(branchEngine, estimatorClones[ip], traceClones[ip], DriftModifier);
        Movers[ip]->initWalkersForPbyP(W.begin() + wPerNode[ip], W.begin() + wPerNode[ip + 1]);
      }
      else
      {
        if (KillNodeCrossing)
          Movers[ip] = new DMCUpdateAllWithKill(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip]);
        else
          Movers[ip] = new DMCUpdateAllWithRejection(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip]);
        Movers[ip]->put(qmcNode);
        Movers[ip]->resetRun(branchEngine, estimatorClones[ip], traceClones[ip], DriftModifier);
        Movers[ip]->initWalkers(W.begin() + wPerNode[ip], W.begin() + wPerNode[ip + 1]);
      }
    }
  }
#if !defined(REMOVE_TRACEMANAGER)
  else
  {
#pragma omp parallel for
    for (int ip = 0; ip < NumThreads; ++ip)
    {
      traceClones[ip]->transfer_state_from(*Traces);
    }
  }
#endif
  branchEngine->checkParameters(W);
  int mxage = mover_MaxAge;
  if (fixW)
  {
    if (BranchInterval < 0)
      BranchInterval = 1;
    mxage = (mover_MaxAge < 0) ? 0 : mover_MaxAge;
    for (int ip = 0; ip < Movers.size(); ++ip)
      Movers[ip]->MaxAge = mxage;
  }
  else
  {
    if (BranchInterval < 0)
      BranchInterval = 1;
    int miage = (qmc_driver_mode[QMC_UPDATE_MODE]) ? 1 : 5;
    mxage     = (mover_MaxAge < 0) ? miage : mover_MaxAge;
    for (int i = 0; i < NumThreads; ++i)
      Movers[i]->MaxAge = mxage;
  }
  {
    std::ostringstream o;
    if (fixW)
      o << "  Fixed population using reconfiguration method\n";
    else
      o << "  Fluctuating population\n";
    o << "  Persistent walkers are killed after " << mxage << " MC sweeps\n";
    o << "  BranchInterval = " << BranchInterval << "\n";
    o << "  Steps per block = " << nSteps << "\n";
    o << "  Number of blocks = " << nBlocks << "\n";
    app_log() << o.str() << std::endl;
  }
  app_log() << "  DMC Engine Initialization = " << init_timer.elapsed() << " secs" << std::endl;
}

bool DMC::run()
{
  LoopTimer dmc_loop;

  bool variablePop = (Reconfiguration == "no");
  resetUpdateEngines();
  //estimator does not need to collect data
  Estimators->setCollectionMode(true);
  Estimators->start(nBlocks);
  for (int ip = 0; ip < NumThreads; ip++)
    Movers[ip]->startRun(nBlocks, false);
#if !defined(REMOVE_TRACEMANAGER)
  Traces->startRun(nBlocks, traceClones);
#endif
  IndexType block        = 0;
  IndexType updatePeriod = (qmc_driver_mode[QMC_UPDATE_MODE]) ? Period4CheckProperties : (nBlocks + 1) * nSteps;
  int sample             = 0;

  RunTimeControl runtimeControl(RunTimeManager, MaxCPUSecs);
  bool enough_time_for_next_iteration = true;

  do // block
  {
    dmc_loop.start();
    Estimators->startBlock(nSteps);
    for (int ip = 0; ip < NumThreads; ip++)
      Movers[ip]->startBlock(nSteps);

    for (IndexType step = 0; step < nSteps; ++step, CurrentStep += BranchInterval)
    {
      //         if(storeConfigs && (CurrentStep%storeConfigs == 0)) {
      //           ForwardWalkingHistory.storeConfigsForForwardWalking(W);
      //           W.resetWalkerParents();
      //         }

#pragma omp parallel
      {
        int ip = omp_get_thread_num();
        Movers[ip]->set_step(sample);
        bool recompute = (step + 1 == nSteps && nBlocksBetweenRecompute && (1 + block) % nBlocksBetweenRecompute == 0 &&
                          qmc_driver_mode[QMC_UPDATE_MODE]);
        wClones[ip]->resetCollectables();
        const size_t nw = W.getActiveWalkers();
#pragma omp for nowait
        for (size_t iw = 0; iw < nw; ++iw)
        {
          Walker_t& thisWalker(*W[iw]);
          Movers[ip]->advanceWalker(thisWalker, recompute);
        }
      }

      //Collectables are weighted but not yet normalized
      if (W.Collectables.size())
      {
        // only when collectable is not empty, need to generalize for W != wClones[0]
        for (int ip = 1; ip < NumThreads; ++ip)
          W.Collectables += wClones[ip]->Collectables;
      }
      branchEngine->branch(CurrentStep, W);
      //         if(storeConfigs && (CurrentStep%storeConfigs == 0)) {
      //           ForwardWalkingHistory.storeConfigsForForwardWalking(W);
      //           W.resetWalkerParents();
      //         }
      if (variablePop)
        FairDivideLow(W.getActiveWalkers(), NumThreads, wPerNode);
      sample++;
    }
    //       branchEngine->debugFWconfig();
    Estimators->stopBlock(acceptRatio());
#if !defined(REMOVE_TRACEMANAGER)
    Traces->write_buffers(traceClones, block);
#endif
    block++;
    if (DumpConfig && block % Period4CheckPoint == 0)
    {
#ifndef USE_FAKE_RNG
      for (int ip = 0; ip < NumThreads; ip++)
        *(RandomNumberControl::Children[ip]) = *(Rng[ip]);
#endif
    }
    recordBlock(block);
    dmc_loop.stop();
    enough_time_for_next_iteration = runtimeControl.enough_time_for_next_iteration(dmc_loop);
    // Rank 0 decides whether the time limit was reached
    myComm->bcast(enough_time_for_next_iteration);

    if (!enough_time_for_next_iteration)
    {
      app_log() << runtimeControl.time_limit_message("DMC", block);
    }
  } while (block < nBlocks && enough_time_for_next_iteration);


  //for(int ip=0; ip<NumThreads; ip++) Movers[ip]->stopRun();
#ifndef USE_FAKE_RNG
  for (int ip = 0; ip < NumThreads; ip++)
    *(RandomNumberControl::Children[ip]) = *(Rng[ip]);
#endif
  Estimators->stop();
  for (int ip = 0; ip < NumThreads; ++ip)
    Movers[ip]->stopRun2();
#if !defined(REMOVE_TRACEMANAGER)
  Traces->stopRun();
#endif
  return finalize(nBlocks);
}


bool DMC::put(xmlNodePtr q)
{
  BranchInterval = -1;
  ParameterSet p;
  p.add(BranchInterval, "branchInterval", "string");
  p.add(BranchInterval, "branchinterval", "string");
  p.add(BranchInterval, "substeps", "int");
  p.add(BranchInterval, "subSteps", "int");
  p.add(BranchInterval, "sub_steps", "int");
  p.put(q);

  //app_log() << "\n DMC::put qmc_counter=" << qmc_common.qmc_counter << "  my_counter=" << MyCounter<< std::endl;
  //app_log() << "  timestep       = " << Tau << std::endl;
  //app_log() << "  blocks         = " << nBlocks << std::endl;
  //app_log() << "  steps          = " << nSteps << std::endl;
  //app_log() << "  current        = " << CurrentStep << std::endl;
  //app_log() << "  walkers/mpi    = " << W.getActiveWalkers() << std::endl << std::endl;
  //app_log().flush();
  return true;
}
} // namespace qmcplusplus
