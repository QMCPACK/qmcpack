//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/DMC/DMC.h"
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/DMC/DMCNonLocalUpdate.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Utilities/Timer.h"

namespace qmcplusplus
{

/// Constructor.
DMC::DMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h), KillNodeCrossing(0), Reconfiguration("no"), Mover(0), BranchInterval(-1)
{
  RootName = "dummy";
  QMCType ="DMC";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_MULTIPLE,1);
  m_param.add(KillWalker,"killnode","string");
  m_param.add(BenchMarkRun,"benchmark","string");
  m_param.add(Reconfiguration,"reconfiguration","string");
  m_param.add(BranchInterval,"branchInterval","int");
  m_param.add(NonLocalMove,"nonlocalmove","string");
  m_param.add(NonLocalMove,"nonlocalmoves","string");
}

bool DMC::run()
{
  resetUpdateEngine();
  //collect is disabled since it is handled by branchEngine->WalkerController
  Estimators->setCollectionMode(true);
  Mover->startRun(nBlocks,true);
  IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
  Timer myclock;
  IndexType block = 0;
  CurrentStep = 0;
  do // block
  {
    Mover->startBlock(nSteps);
    for(IndexType step=0; step< nSteps; step++, CurrentStep+=BranchInterval)
    {
      if(storeConfigs && (CurrentStep%storeConfigs == 0))
      {
        ForwardWalkingHistory.storeConfigsForForwardWalking(W);
        W.resetWalkerParents();
      }
      for(IndexType interval=0; interval<BranchInterval-1; ++interval)
        Mover->advanceWalkers(W.begin(),W.end(), false);
      W.resetCollectables();
      Mover->advanceWalkers(W.begin(),W.end(), false);
      Mover->setMultiplicity(W.begin(),W.end());
      branchEngine->branch(CurrentStep,W);
    }
    block++;
    Mover->stopBlock();
    recordBlock(block);
    if(QMCDriverMode[QMC_UPDATE_MODE] && CurrentStep%updatePeriod == 0)
      Mover->updateWalkers(W.begin(), W.end());
  }
  while(block<nBlocks &&  myclock.elapsed()<MaxCPUSecs);
  Mover->stopRun();
  return finalize(block);
}

void DMC::resetUpdateEngine()
{
  bool fixW=((Reconfiguration == "yes")||(Reconfiguration == "pure"));
  if(Mover==0) //disable switching update modes for DMC in a run
  {
    //load walkers if they were saved
    W.loadEnsemble();
    branchEngine->initWalkerController(W,fixW);
    if(QMCDriverMode[QMC_UPDATE_MODE])
    {
      if(NonLocalMove == "yes")
      {
        app_log() << "  Non-local update is used." << std::endl;
        DMCNonLocalUpdatePbyP* nlocMover= new DMCNonLocalUpdatePbyP(W,Psi,H,Random);
        nlocMover->put(qmcNode);
        Mover=nlocMover;
      }
      else
      {
        Mover= new DMCUpdatePbyPWithRejection(W,Psi,H,Random);
      }
      Mover->resetRun(branchEngine,Estimators);
      Mover->initWalkersForPbyP(W.begin(),W.end());
    }
    else
    {
      if(NonLocalMove == "yes")
      {
        app_log() << "  Non-local update is used." << std::endl;
        DMCNonLocalUpdate* nlocMover= new DMCNonLocalUpdate(W,Psi,H,Random);
        nlocMover->put(qmcNode);
        Mover=nlocMover;
      }
      else
      {
        if(KillNodeCrossing)
        {
          Mover = new DMCUpdateAllWithKill(W,Psi,H,Random);
        }
        else
        {
          Mover = new DMCUpdateAllWithRejection(W,Psi,H,Random);
        }
      }
      Mover->resetRun(branchEngine,Estimators);
      Mover->initWalkers(W.begin(),W.end());
    }
    branchEngine->checkParameters(W);
  }
  else
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
      Mover->updateWalkers(W.begin(),W.end());
  }
  if(fixW)
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
      app_log() << "  DMC PbyP Update with reconfigurations" << std::endl;
    else
      app_log() << "  DMC walker Update with reconfigurations" << std::endl;
    Mover->MaxAge=0;
    if(BranchInterval<0)
    {
      BranchInterval=nSteps;
    }
    //nSteps=1;
  }
  else
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
    {
      app_log() << "  DMC PbyP Update with a fluctuating population" << std::endl;
      Mover->MaxAge=1;
    }
    else
    {
      app_log() << "  DMC walker Update with a fluctuating population" << std::endl;
      Mover->MaxAge=3;
    }
    if(BranchInterval<0)
      BranchInterval=1;
  }
  app_log() << "  BranchInterval = " << BranchInterval << std::endl;
  app_log() << "  Steps per block = " << nSteps << std::endl;
  app_log() << "  Number of blocks = " << nBlocks << std::endl;
}

bool DMC::put(xmlNodePtr q)
{
  return true;
}

}

