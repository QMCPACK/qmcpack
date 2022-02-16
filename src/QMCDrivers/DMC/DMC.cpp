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


#include "DMC.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/DMC/DMCUpdatePbyPL2.h"
#include "QMCDrivers/DMC/SODMCUpdatePbyP.h"
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Concurrency/OpenMP.h"
#include "Utilities/Timer.h"
#include "Utilities/RunTimeManager.h"
#include "RandomNumberControl.h"
#include "Utilities/ProgressReportEngine.h"
#include "Utilities/qmc_common.h"
#include "Utilities/FairDivide.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
using TraceManager = int;
#endif

namespace qmcplusplus
{
/// Constructor.
DMC::DMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, Communicate* comm, bool enable_profiling)
    : QMCDriver(w, psi, h, comm, "DMC", enable_profiling),
      KillNodeCrossing(0),
      BranchInterval(-1),
      L2("no"),
      Reconfiguration("no"),
      mover_MaxAge(-1)
{
  RootName = "dmc";
  qmc_driver_mode.set(QMC_UPDATE_MODE, 1);
  m_param.add(KillWalker, "killnode");
  m_param.add(Reconfiguration, "reconfiguration");
  //m_param.add(BranchInterval,"branchInterval");
  m_param.add(NonLocalMove, "nonlocalmove");
  m_param.add(NonLocalMove, "nonlocalmoves");
  m_param.add(mover_MaxAge, "MaxAge");
  m_param.add(L2, "L2_diffusion");
}

void DMC::resetUpdateEngines()
{
  ReportEngine PRE("DMC", "resetUpdateEngines");
  bool fixW = (Reconfiguration == "runwhileincorrect");
  if (Reconfiguration != "no" && Reconfiguration != "runwhileincorrect")
    throw std::runtime_error("Reconfiguration is currently broken and gives incorrect results. Use dynamic "
                             "population control by setting reconfiguration=\"no\" or removing the reconfiguration "
                             "option from the DMC input section. If accessing the broken reconfiguration code path "
                             "is still desired, set reconfiguration to \"runwhileincorrect\" instead of \"yes\".");
  makeClones(W, Psi, H);
  Timer init_timer;
  bool spinor = false;
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
    Rng.resize(NumThreads);
    estimatorClones.resize(NumThreads, 0);
    traceClones.resize(NumThreads, 0);
    FairDivideLow(W.getActiveWalkers(), NumThreads, wPerRank);

    {
      //log file
      std::ostringstream o;
      o << "  Initial partition of walkers on a node: ";
      copy(wPerRank.begin(), wPerRank.end(), std::ostream_iterator<int>(o, " "));
      o << "\n";
      if (qmc_driver_mode[QMC_UPDATE_MODE])
      {
        o << "  Updates by particle-by-particle moves";
        if (L2 == "yes")
          app_log() << "Using DMCUpdatePbyPL2" << std::endl;
        else
          app_log() << "Using DMCUpdatePbyPWithRejectionFast" << std::endl;
      }
      else
        o << "  Updates by walker moves";
      // Appears to be set in constructor reported here and used nowhere
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
      Rng[ip] = std::make_unique<FakeRandom>();
#else
      Rng[ip] = std::make_unique<RandomGenerator>(*RandomNumberControl::Children[ip]);
      hClones[ip]->setRandomGenerator(Rng[ip].get());
#endif
      if (W.isSpinor())
      {
        spinor = true;
        if (qmc_driver_mode[QMC_UPDATE_MODE])
        {
          Movers[ip] = new SODMCUpdatePbyPWithRejectionFast(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip]);
          Movers[ip]->setSpinMass(SpinMass);
          Movers[ip]->put(qmcNode);
          Movers[ip]->resetRun(branchEngine.get(), estimatorClones[ip], traceClones[ip], DriftModifier);
          Movers[ip]->initWalkersForPbyP(W.begin() + wPerRank[ip], W.begin() + wPerRank[ip + 1]);
        }
        else
        {
          APP_ABORT("SODMC Driver Mode must be PbyP\n");
        }
      }
      else
      {
        if (qmc_driver_mode[QMC_UPDATE_MODE])
        {
          if (L2 == "yes")
            Movers[ip] = new DMCUpdatePbyPL2(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip]);
          else
            Movers[ip] = new DMCUpdatePbyPWithRejectionFast(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip]);

          Movers[ip]->put(qmcNode);
          Movers[ip]->resetRun(branchEngine.get(), estimatorClones[ip], traceClones[ip], DriftModifier);
          Movers[ip]->initWalkersForPbyP(W.begin() + wPerRank[ip], W.begin() + wPerRank[ip + 1]);
        }
        else
        {
          if (KillNodeCrossing)
            Movers[ip] = new DMCUpdateAllWithKill(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip]);
          else
            Movers[ip] = new DMCUpdateAllWithRejection(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip]);
          Movers[ip]->put(qmcNode);
          Movers[ip]->resetRun(branchEngine.get(), estimatorClones[ip], traceClones[ip], DriftModifier);
          Movers[ip]->initWalkers(W.begin() + wPerRank[ip], W.begin() + wPerRank[ip + 1]);
        }
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

  if (spinor)
    app_log() << "   Spins treated as dynamic variable with SpinMass: " << SpinMass << std::endl;

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
  LoopTimer<> dmc_loop;

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

  RunTimeControl<> runtimeControl(run_time_manager, MaxCPUSecs, myComm->getName(), myComm->rank() == 0);

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
        FairDivideLow(W.getActiveWalkers(), NumThreads, wPerRank);
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
        *RandomNumberControl::Children[ip] = *Rng[ip];
#endif
    }
    recordBlock(block);
    dmc_loop.stop();

    bool stop_requested = false;
    // Rank 0 decides whether the time limit was reached
    if (!myComm->rank())
      stop_requested = runtimeControl.checkStop(dmc_loop);
    myComm->bcast(stop_requested);

    if (stop_requested)
    {
      if (!myComm->rank())
        app_log() << runtimeControl.generateStopMessage("DMC", block - 1);
      run_time_manager.markStop();
      break;
    }

  } while (block < nBlocks);

#ifndef USE_FAKE_RNG
  for (int ip = 0; ip < NumThreads; ip++)
    *RandomNumberControl::Children[ip] = *Rng[ip];
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
  p.add(BranchInterval, "branchInterval");
  p.add(BranchInterval, "branchinterval");
  p.add(BranchInterval, "substeps");
  p.add(BranchInterval, "subSteps");
  p.add(BranchInterval, "sub_steps");
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
