//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/DMC/DMCPbyP.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/DMC/DMCNonLocalUpdate.h"

namespace qmcplusplus { 

  /// Constructor.
  DMCPbyP::DMCPbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h),
    KillNodeCrossing(0),
    PopIndex(-1), EtrialIndex(-1),BranchInterval(-1), NonLocalMoveIndex(-1),
    BranchInfo("default"), KillWalker("no"), Reconfiguration("no"), NonLocalMove("no"), 
    Mover(0){ 
    RootName = "dmc";
    QMCType ="DMCPbyP";

    //to prevent initialization
    QMCDriverMode.set(QMC_MULTIPLE,1);
    QMCDriverMode.set(QMC_UPDATE_MODE,1);

    m_param.add(KillWalker,"killnode","string");
    m_param.add(Reconfiguration,"reconfiguration","string");
    m_param.add(BranchInterval,"branch_interval","int");
    m_param.add(BranchInterval,"branchInterval","int");
    m_param.add(NonLocalMove,"nonlocalmoves","string");
  }
  
  /// destructor
  DMCPbyP::~DMCPbyP() { 
    if(Mover) delete Mover;
  }

  bool DMCPbyP::run() { 

    bool fixW = (Reconfiguration == "yes");
    bool killNC = (KillWalker == "yes");

    branchEngine->initWalkerController(Tau,fixW);

    if(Mover ==0) {
      if(NonLocalMove == "yes") {
        app_log() << "  Non-local update is used." << endl;
        DMCNonLocalUpdatePbyP* nlocMover= new DMCNonLocalUpdatePbyP(W,Psi,H,Random);
        nlocMover->put(qmcNode);
        Mover=nlocMover;
        NonLocalMoveIndex=Estimators->addColumn("NonLocalMove");
      } else {
        if(killNC) {
          app_log() << "  Kill a walker, if a node crossing is detected." << endl;
          Mover = new DMCUpdatePbyPWithKill(W,Psi,H,Random);
        } else {
          app_log() << "  Reject a move when the node crossing is detected." << endl;
          Mover = new DMCUpdatePbyPWithRejection(W,Psi,H,Random);
        }
      }
    }

    //set the collection mode for the estimator
    Estimators->setCollectionMode(branchEngine->SwapMode);

    PopIndex = Estimators->addColumn("Population");
    EtrialIndex = Estimators->addColumn("Etrial");
    Estimators->reportHeader(AppendRun);
    Estimators->reset();

    Mover->resetRun(branchEngine);
    Mover->initWalkers(W.begin(),W.end());

    if(fixW)  {
      if(BranchInterval<0) {
        BranchInterval=nSteps;
        nSteps=1;
      }
      app_log() << "  DMC PbyP update with reconfigurations" << endl;
      app_log() << "    BranchInterval=" << BranchInterval << endl;
      app_log() << "    Steps         =" << nSteps << endl;
      app_log() << "    Blocks        =" << nBlocks << endl;
      dmcWithReconfiguration();
    } else {
      app_log() << "  DMC PbyP update with a fluctuating population" << endl;
      dmcWithBranching();
    }

    Estimators->finalize();

    return true;
  }

  void DMCPbyP::dmcWithReconfiguration() {
    //MaxAge is set to 0 not to evaluate branching factor
    bool checkNonLocalMove=(NonLocalMoveIndex>0);
    Mover->MaxAge=0;
    IndexType block = 0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;
    RealType Eest=branchEngine->E_T;
    do {
      IndexType step = 0;
      Mover->startBlock();
      Estimators->startBlock();
      do {
        int interval=0; 
        do {
          Mover->advanceWalkers(W.begin(), W.end());
          ++interval;
          ++CurrentStep; 
        } while(interval<BranchInterval);

        Estimators->accumulate(W);
        int nwKept= branchEngine->branch(CurrentStep,W);
        Estimators->setColumn(PopIndex,nwKept);
        Eest = branchEngine->CollectAndUpdate(W.getActiveWalkers(),Eest);

        ++step; 
      } while(step<nSteps);

      if(checkNonLocalMove) {
        Estimators->setColumn(NonLocalMoveIndex, 
            static_cast<RealType>(Mover->NonLocalMoveAccepted)/static_cast<RealType>(W.getActiveWalkers()*BranchInterval));
      }

      Estimators->stopBlock(Mover->acceptRatio());

      nAcceptTot += Mover->nAccept;
      nRejectTot += Mover->nReject;
      
      Estimators->setColumn(EtrialIndex,branchEngine->E_T);
      block++;
      recordBlock(block);

      updateWalkers();
    } while(block<nBlocks);
  }

  void DMCPbyP::dmcWithBranching() {

    bool checkNonLocalMove=(NonLocalMoveIndex>0);
    Mover->MaxAge=1;
    IndexType block = 0;
    RealType Eest = branchEngine->E_T;
    nAcceptTot = 0;
    nRejectTot = 0;

    do {
      IndexType step = 0;
      IndexType pop_acc=0; 

      Mover->startBlock();
      Estimators->startBlock();
      do {

        IndexType interval = 0;
        do {
          Mover->advanceWalkers(W.begin(),W.end());
          ++interval;
          ++CurrentStep;
        } while(interval<BranchInterval);

        ++step; 
        Mover->setMultiplicity(W.begin(),W.end());

        //Will merge
        Estimators->accumulate(W);
        int cur_pop = branchEngine->branch(CurrentStep,W);
        Eest = branchEngine->CollectAndUpdate(cur_pop,Eest);

        pop_acc += cur_pop;
        if(CurrentStep%100 == 0) updateWalkers();
      } while(step<nSteps);
      
      if(checkNonLocalMove) {
        Estimators->setColumn(NonLocalMoveIndex, 
        static_cast<RealType>(Mover->NonLocalMoveAccepted)/
        static_cast<RealType>(W.getActiveWalkers()*BranchInterval));
      }

      Estimators->stopBlock(Mover->acceptRatio());

      nAcceptTot += Mover->nAccept;
      nRejectTot += Mover->nReject;

      //update estimator
      Estimators->setColumn(PopIndex,static_cast<RealType>(pop_acc)/static_cast<RealType>(nSteps));
      Estimators->setColumn(EtrialIndex,Eest); 
      Eest = Estimators->average(0);
      RealType totmoves=1.0/static_cast<RealType>(step*W.getActiveWalkers());

      block++;
      recordBlock(block);
    } while(block<nBlocks);
  }

  bool DMCPbyP::put(xmlNodePtr q){
    return true;
  }

}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
