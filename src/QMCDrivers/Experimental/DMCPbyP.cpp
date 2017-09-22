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
    
    
#include "QMCDrivers/DMC/DMCPbyP.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/DMC/DMCNonLocalUpdate.h"
#include "Estimators/DMCEnergyEstimator.h"

namespace qmcplusplus
{

/// Constructor.
DMCPbyP::DMCPbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h),
  KillNodeCrossing(0),
  BranchInterval(-1), NonLocalMoveIndex(-1),
  BranchInfo("default"), KillWalker("no"), Reconfiguration("no"), NonLocalMove("no"),
  Mover(0)
{
  RootName = "dmc";
  QMCType ="DMCPbyP";
  //to prevent initialization
  QMCDriverMode.set(QMC_MULTIPLE,1);
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  m_param.add(KillWalker,"killnode","string");
  m_param.add(Reconfiguration,"reconfiguration","string");
  m_param.add(BranchInterval,"branch_interval","int");
  m_param.add(BranchInterval,"branchInterval","int");
  m_param.add(NonLocalMove,"nonlocalmove","string");
  m_param.add(NonLocalMove,"nonlocalmoves","string");
  //create a ScalarEstimator and add DMCEnergyEstimator
  Estimators = new EstimatorManagerBase(H);
  Estimators->add(new DMCEnergyEstimator,"elocal");
}

/// destructor
DMCPbyP::~DMCPbyP()
{
  if(Mover)
    delete Mover;
}

bool DMCPbyP::run()
{
  bool fixW = (Reconfiguration == "yes");
  bool killNC = (KillWalker == "yes");
  branchEngine->initWalkerController(Tau,fixW);
  if(Mover ==0)
  {
    if(NonLocalMove == "yes")
    {
      app_log() << "  Non-local update is used." << std::endl;
      DMCNonLocalUpdatePbyP* nlocMover= new DMCNonLocalUpdatePbyP(W,Psi,H,Random);
      nlocMover->put(qmcNode);
      Mover=nlocMover;
      NonLocalMoveIndex=Estimators->addColumn("NonLocalMove");
    }
    else
    {
      if(killNC)
      {
        app_log() << "  Kill a walker, if a node crossing is detected." << std::endl;
        Mover = new DMCUpdatePbyPWithKill(W,Psi,H,Random);
      }
      else
      {
        app_log() << "  Reject a move when the node crossing is detected." << std::endl;
        Mover = new DMCUpdatePbyPWithRejection(W,Psi,H,Random);
      }
    }
  }
  //set the collection mode for the estimator
  Estimators->setCollectionMode(branchEngine->SwapMode);
  Estimators->reportHeader(AppendRun);
  Estimators->reset();
  Mover->resetRun(branchEngine);
  Mover->initWalkers(W.begin(),W.end());
  if(fixW)
  {
    Mover->MaxAge=0;
    if(BranchInterval<0)
    {
      BranchInterval=nSteps;
      nSteps=1;
    }
    app_log() << "  DMC PbyP update with reconfigurations" << std::endl;
  }
  else
  {
    Mover->MaxAge=1;
    if(BranchInterval<0)
      BranchInterval=1;
    app_log() << "  DMC PbyP update with a fluctuating population" << std::endl;
  }
  app_log() << "    BranchInterval=" << BranchInterval << std::endl;
  app_log() << "    Steps         =" << nSteps << std::endl;
  app_log() << "    Blocks        =" << nBlocks << std::endl;
  nAcceptTot = 0;
  nRejectTot = 0;
  return dmcWithBranching();
}

bool DMCPbyP::dmcWithBranching()
{
  bool checkNonLocalMove=(NonLocalMoveIndex>0);
  //Mover->MaxAge=1;
  IndexType block = 0;
  RealType Eest = branchEngine->E_T;
  do
  {
    IndexType step = 0;
    IndexType pop_acc=0;
    Mover->startBlock();
    Estimators->startBlock();
    Mover->NonLocalMoveAccepted=0;
    do
    {
      //DMC without branching, weights are accumulated
      IndexType interval = 0;
      do
      {
        Mover->advanceWalkers(W.begin(),W.end());
        ++interval;
        ++CurrentStep;
      }
      while(interval<BranchInterval);
      //set the multiplicity with the weights
      Mover->setMultiplicity(W.begin(),W.end());
      Estimators->accumulate(W);
      branchEngine->branch(CurrentStep,W);
      pop_acc += W.getActiveWalkers();
//         if(CurrentStep%100 == 0) updateWalkers();
      if(CurrentStep%Period4CheckProperties  == 0)
        updateWalkers();
      ++step;
    }
    while(step<nSteps);
    if(checkNonLocalMove)
    {
      Estimators->setColumn(NonLocalMoveIndex,
                            static_cast<RealType>(Mover->NonLocalMoveAccepted)/
                            static_cast<RealType>(W.getActiveWalkers()*BranchInterval));
    }
    Estimators->stopBlock(Mover->acceptRatio());
    nAcceptTot += Mover->nAccept;
    nRejectTot += Mover->nReject;
    Eest = Estimators->average(0);
    RealType totmoves=1.0/static_cast<RealType>(step*W.getActiveWalkers());
    block++;
    recordBlock(block);
  }
  while(block<nBlocks);
  return finalize(block);
}

bool DMCPbyP::dmcWithReconfiguration()
{
  //MaxAge is set to 0 not to evaluate branching factor
  bool checkNonLocalMove=(NonLocalMoveIndex>0);
  IndexType block = 0;
  RealType Eest=branchEngine->E_T;
  do
  {
    IndexType step = 0;
    Mover->startBlock();
    Estimators->startBlock();
    do
    {
      int interval=0;
      Mover->NonLocalMoveAccepted=0;
      do
      {
        Mover->advanceWalkers(W.begin(), W.end());
        ++interval;
        ++CurrentStep;
      }
      while(interval<BranchInterval);
      Estimators->accumulate(W);
      branchEngine->branch(CurrentStep,W);
      ++step;
    }
    while(step<nSteps);
    if(checkNonLocalMove)
    {
      Estimators->setColumn(NonLocalMoveIndex,
                            static_cast<RealType>(Mover->NonLocalMoveAccepted)/static_cast<RealType>(W.getActiveWalkers()*BranchInterval));
    }
    Estimators->stopBlock(Mover->acceptRatio());
    nAcceptTot += Mover->nAccept;
    nRejectTot += Mover->nReject;
    block++;
    recordBlock(block);
    updateWalkers();
  }
  while(block<nBlocks);
  return finalize(block);
}

bool DMCPbyP::put(xmlNodePtr q)
{
  return true;
}

}

