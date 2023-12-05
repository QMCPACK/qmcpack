//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "RMC.h"
#include "QMCDrivers/RMC/RMCUpdatePbyP.h"
#include "QMCDrivers/RMC/RMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "RandomNumberControl.h"
#include "Concurrency/OpenMP.h"
#include "Message/CommOperators.h"
#include "Particle/Reptile.h"
#include "Utilities/FairDivide.h"
#include "Utilities/RunTimeManager.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
using TraceManager = int;
#endif


namespace qmcplusplus
{
/// Constructor.
RMC::RMC(const ProjectData& project_data,
         MCWalkerConfiguration& w,
         TrialWaveFunction& psi,
         QMCHamiltonian& h,
         Communicate* comm)
    : QMCDriver(project_data, w, psi, h, comm, "RMC"),
      prestepsVMC(-1),
      rescaleDrift("no"),
      beta(-1),
      beads(-1),
      fromScratch(true)
{
  RootName = "rmc";
  qmc_driver_mode.set(QMC_UPDATE_MODE, 1);
  qmc_driver_mode.set(QMC_WARMUP, 0);
  m_param.add(rescaleDrift, "drift");
  m_param.add(beta, "beta");
  m_param.add(beads, "beads");
  m_param.add(resizeReptile, "resize");
  m_param.add(prestepsVMC, "vmcpresteps");

  Action.resize(3);
  Action[0] = w.addProperty("ActionBackward");
  Action[1] = w.addProperty("ActionForward");
  Action[2] = w.addProperty("ActionLocal");
  TransProb.resize(2);
  TransProb[0] = w.addProperty("TransProbBackward");
  TransProb[1] = w.addProperty("TransProbForward");
}

bool RMC::run()
{
  resetRun();
  //start the main estimator
  Estimators->start(nBlocks);
  for (int ip = 0; ip < NumThreads; ++ip)
    Movers[ip]->startRun(nBlocks, false);
#if !defined(REMOVE_TRACEMANAGER)
  Traces->startRun(nBlocks, traceClones);
#endif
  const bool has_collectables = W.Collectables.size();

  LoopTimer<> rmc_loop;
  RunTimeControl<> runtimeControl(run_time_manager, MaxCPUSecs, myComm->getName(), myComm->rank() == 0);
  for (int block = 0; block < nBlocks; ++block)
  {
    rmc_loop.start();
#pragma omp parallel
    {
      int ip                 = omp_get_thread_num();
      IndexType updatePeriod = (qmc_driver_mode[QMC_UPDATE_MODE]) ? Period4CheckProperties : 0;
      //assign the iterators and resuse them
      MCWalkerConfiguration::iterator wit(W.begin() + wPerRank[ip]), wit_end(W.begin() + wPerRank[ip + 1]);
      Movers[ip]->startBlock(nSteps);
      int now_loc = CurrentStep;

      RealType cnorm = 1.0; //This is because there is only one reptile per walkerset.

      for (int step = 0; step < nSteps; ++step)
      {
        //collectables are reset, it is accumulated while advancing walkers
        wClones[ip]->resetCollectables();
        Movers[ip]->advanceWalkers(wit, wit_end, false);
        if (has_collectables)
          wClones[ip]->Collectables *= cnorm;
        Movers[ip]->accumulate(wit, wit_end);

        ++now_loc;
        if (Period4WalkerDump && now_loc % myPeriod4WalkerDump == 0)
          wClones[ip]->saveEnsemble(wit, wit_end);

        branchEngine->collect(
            CurrentStep,
            W); //Ray Clay:  For now, collects and syncs based on first reptile.  Need a better way to do this.
      }
      Movers[ip]->stopBlock(false);
    } //end-of-parallel for
    CurrentStep += nSteps;
    Estimators->stopBlock(estimatorClones);
    recordBlock(block);
    rmc_loop.stop();

    bool stop_requested = false;
    // Rank 0 decides whether the time limit was reached
    if (!myComm->rank())
      stop_requested = runtimeControl.checkStop(rmc_loop);
    myComm->bcast(stop_requested);

    if (stop_requested)
    {
      if (!myComm->rank())
        app_log() << runtimeControl.generateStopMessage("RMC", block);
      run_time_manager.markStop();
      break;
    }
  } //block
  Estimators->stop(estimatorClones);
  //copy back the random states
  for (int ip = 0; ip < NumThreads; ++ip)
    RandomNumberControl::Children[ip] = Rng[ip]->makeClone();
  //return nbeads and stuff to its original unset state;
  resetVars();
  return finalize(nBlocks);
}

void RMC::resetRun()
{
  m_param.put(qmcNode);
  //For now, assume that nReptiles=NumThreads;
  nReptiles = NumThreads;

  if (beads < 1)
    beads = beta / Tau;
  else
    beta = beads * Tau;

  app_log() << "Projection time:  " << beta << " Ha^-1" << std::endl;
  //Calculate the number of VMC presteps if not given:
  if (prestepsVMC == -1 && fromScratch == true)
    prestepsVMC = beads + 2;
  //Check to see if the MCWalkerConfiguration is in a state suitable for reptation
  if (!W.ReptileList.empty())
  {
    fromScratch = false;

    app_log() << "Previous RMC reptiles detected...\n";
    if (Tau == W.ReptileList[0]->getTau() && beads == W.ReptileList[0]->size())
      app_log() << "  Using current reptiles\n"; //do nothing
    else                                         //we need to extrapolate off of the current reptile set.
    {
      //pull the reptile configurations out
      app_log() << "  Previous Tau/Beta:  " << W.ReptileList[0]->getTau() << "/"
                << W.ReptileList[0]->getTau() * W.ReptileList[0]->size() << std::endl;
      app_log() << "  New      Tau/Beta: " << Tau << "/" << beta << std::endl;
      app_log() << "    Linear interpolation to get new reptile.\n";
      std::vector<ReptileConfig_t> repSamps(0);
      for (IndexType sampid = 0; sampid < W.ReptileList.size() && sampid < nReptiles; sampid++)
        repSamps.push_back(W.ReptileList[sampid]->getReptileSlicePositions(Tau, beta));

      //In the event of a disparity in the number of requested reptiles and the ones received....  just copy
      //Copies cyclically.  First iteration copies the first entry, second the second, and so on.  So we don't replicate just one config.
      for (IndexType copyid = 0; repSamps.size() < nReptiles; copyid++)
        repSamps.push_back(repSamps[copyid]);


      resetReptiles(repSamps, Tau);
    }
  }

  //Previous run was nothing, VMC, or DMC.  No reptiles--so we initialize based on whatever is there.
  else
  {
    //Initialize on whatever walkers are in MCWalkerConfiguration.
    app_log() << "Using walkers from previous non-RMC run.\n";
    std::vector<ParticlePos> wSamps(0);
    MCWalkerConfiguration::iterator wit(W.begin()), wend(W.end());
    for (IndexType sampid = 0; wit != wend && sampid < nReptiles; wit++)
      wSamps.push_back((**wit).R);

    for (IndexType copyid = 0; wSamps.size() < nReptiles; copyid++)
      wSamps.push_back(wSamps[copyid]);
    resetReptiles(wSamps, beads, Tau);
  }

  //Now that we know if we're starting from scratch... decide whether to force VMC warmup.
  if (prestepsVMC == -1 && fromScratch == true)
    prestepsVMC = beads + 2;
  makeClones(W, Psi, H);
  myPeriod4WalkerDump = (Period4WalkerDump > 0) ? Period4WalkerDump : (nBlocks + 1) * nSteps;

  if (Movers.empty())
  {
    Movers.resize(NumThreads, nullptr);
    estimatorClones.resize(NumThreads, nullptr);
    traceClones.resize(NumThreads, nullptr);
    Rng.resize(NumThreads);
    branchEngine->initReptile(W);

    // hdf_archive::hdf_archive() is not thread-safe
    for (int ip = 0; ip < NumThreads; ++ip)
      estimatorClones[ip] = new EstimatorManagerBase(*Estimators);

#pragma omp parallel for
    for (int ip = 0; ip < NumThreads; ++ip)
    {
      std::ostringstream os;
      estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
      estimatorClones[ip]->setCollectionMode(false);
      Rng[ip] = RandomNumberControl::Children[ip]->makeClone();
#if !defined(REMOVE_TRACEMANAGER)
      traceClones[ip] = Traces->makeClone();
#endif
      hClones[ip]->setRandomGenerator(Rng[ip].get());
      if (qmc_driver_mode[QMC_UPDATE_MODE])
      {
        os << "  PbyP moves with drift, using RMCUpdatePbyPWithDriftFast" << std::endl;
        Movers[ip] =
            new RMCUpdatePbyPWithDrift(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip], Action, TransProb);
      }
      else
      {
        os << "  walker moves with drift, using RMCUpdateAllWithDriftFast" << std::endl;
        Movers[ip] = new RMCUpdateAllWithDrift(*wClones[ip], *psiClones[ip], *hClones[ip], *Rng[ip], Action, TransProb);
      }
      Movers[ip]->nSubSteps = nSubSteps;
      if (ip == 0)
        app_log() << os.str() << std::endl;
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
  app_log().flush();
#pragma omp parallel for
  for (int ip = 0; ip < NumThreads; ++ip)
  {
    Movers[ip]->put(qmcNode);
    Movers[ip]->resetRun(branchEngine.get(), estimatorClones[ip], traceClones[ip], DriftModifier);

    wClones[ip]->reptile    = W.ReptileList[ip].get();
    wClones[ip]->activeBead = 0;
    wClones[ip]->direction  = +1;

    if (qmc_driver_mode[QMC_UPDATE_MODE])
    {
      // app_log () << ip << " initWalkers for pbyp...\n";
      Movers[ip]->initWalkersForPbyP(W.ReptileList[ip]->repstart, W.ReptileList[ip]->repend);
    }
    else
    {
      Movers[ip]->initWalkers(W.begin() + wPerRank[ip], W.begin() + wPerRank[ip + 1]);
    }

    //this will "unroll" the reptile according to forced VMC steps (no bounce).  See beginning of function for logic of setting prestepVMC.
    for (IndexType prestep = 0; prestep < prestepsVMC; prestep++)
    {
      Movers[ip]->advanceWalkers(W.begin(), W.begin(), true);
    }

    //set up initial action and transprob.
    MCWalkerConfiguration::iterator wit(W.begin() + wPerRank[ip]), wit_end(W.begin() + wPerRank[ip + 1]);
  }


  app_log() << "Finished " << prestepsVMC << " VMC presteps\n";
  branchEngine->checkParameters(W);

#pragma omp parallel for
  for (int ip = 0; ip < NumThreads; ++ip)
  {
    for (int prestep = 0; prestep < nWarmupSteps; ++prestep)
    {
      Movers[ip]->advanceWalkers(W.begin() + wPerRank[ip], W.begin() + wPerRank[ip + 1], false);
      branchEngine->collect(CurrentStep, W);
    }
  }

  fromScratch = false;
}

bool RMC::put(xmlNodePtr q)
{
  m_param.put(q);
  return true;
}

//This will resize the MCWalkerConfiguration and initialize the ReptileList.  Does not care for previous runs.
void RMC::resetReptiles(int nReptiles_in, int nbeads_in, RealType tau)
{
  W.ReptileList.clear();
  // Maybe we should be more vigorous in cleaning the MCWC WalkerList?
  std::vector<int> repWalkerSlice;
  int nwtot = nbeads_in * nReptiles_in;
  FairDivideLow(nwtot, nReptiles_in, repWalkerSlice);
  if (W.getActiveWalkers() - nwtot != 0)
    addWalkers(nwtot - W.getActiveWalkers());

  for (int i = 0; i < nReptiles_in; i++)
  {
    W.ReptileList.push_back(
        std::make_unique<Reptile>(W, W.begin() + repWalkerSlice[i], W.begin() + repWalkerSlice[i + 1]));
    W.ReptileList[i]->setTau(tau);
  }
}
//This will resize the MCWalkerConfiguration and initialize Reptile list.  It will then reinitialize the MCWC with a list of Reptile coordinates
void RMC::resetReptiles(std::vector<ReptileConfig_t>& reptile_samps, RealType tau)
{
  if (reptile_samps.empty())
  {
    APP_ABORT("RMC::resetReptiles(std::vector< ReptileConfig_t > reptile_samps):  No samples!\n");
  }
  else
  {
    IndexType nReptiles_in = reptile_samps.size();
    IndexType nBeads_in    = reptile_samps[0].size();
    resetReptiles(nReptiles_in, nBeads_in, tau);

    for (IndexType i = 0; i < W.ReptileList.size(); i++)
    {
      W.ReptileList[i]->setReptileSlicePositions(reptile_samps[i]);
    }
  }
}
//For # of walker samples, create that many reptiles with nbeads each.  Initialize each reptile to have the value of the walker "seed".
void RMC::resetReptiles(std::vector<ParticlePos>& walker_samps, int nBeads_in, RealType tau)
{
  if (walker_samps.empty())
  {
    APP_ABORT("RMC::resetReptiles(std::vector< ParticlePos > walker_samps):  No samples!\n");
  }
  else
  {
    IndexType nReptiles_in = walker_samps.size();
    resetReptiles(nReptiles_in, nBeads_in, tau);

    for (IndexType i = 0; i < W.ReptileList.size(); i++)
    {
      W.ReptileList[i]->setReptileSlicePositions(walker_samps[i]);
    }
  }
}

}; // namespace qmcplusplus
