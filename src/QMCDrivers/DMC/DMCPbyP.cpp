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

namespace qmcplusplus { 

  /// Constructor.
  DMCPbyP::DMCPbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h),
    KillNodeCrossing(0),
    PopIndex(-1), EtrialIndex(-1),
    BranchInfo("default"), KillWalker("no"), Reconfiguration("no"), Mover(0){ 
    RootName = "dmc";
    QMCType ="dmc";

    //to prevent initialization
    QMCDriverMode.set(QMC_MULTIPLE,1);
    QMCDriverMode.set(QMC_UPDATE_MODE,1);

    m_param.add(KillWalker,"killnode","string");
    m_param.add(Reconfiguration,"reconfiguration","string");
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
      //if(killNC) {
      //  app_log() << "Walkers will be killed if a node crossing is detected." << endl;
      //} else {
      //  app_log() << "Walkers will be kept even if a node crossing is detected." << endl;
      //}
      app_log() << "  Reject a move when the node crossing is detected." << endl;
      Mover = new DMCUpdatePbyPWithRejection(W,Psi,H,Random);
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
      app_log() << "  DMC Particle-By-Particle Update with reconfigurations" << endl;
      dmcWithReconfiguration();
    } else {
      app_log() << "  DMC Particle-By-Particle Update with a fluctuating population" << endl;
      dmcWithBranching();
    }

    Estimators->finalize();

    return true;
  }

  void DMCPbyP::dmcWithReconfiguration() {
    IndexType block = 0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;
    do {
      IndexType step = 0;
      Mover->startBlock();
      Estimators->startBlock();
      do {
        Mover->advanceWalkers(W.begin(), W.end());
        MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
        while(it != it_end) {
          branchEngine->accumulate((*it)->Properties(LOCALENERGY),1.0);
          ++it;
        }
        step++; CurrentStep++;
      } while(step<nSteps);

      Estimators->accumulate(W);
      int nwKept= branchEngine->branch(CurrentStep,W);

      Estimators->stopBlock(Mover->acceptRatio());

      nAcceptTot += Mover->nAccept;
      nRejectTot += Mover->nReject;
      
      Estimators->setColumn(PopIndex,nwKept);
      Estimators->setColumn(EtrialIndex,branchEngine->E_T);

      block++;
      recordBlock(block);

      updateWalkers();
    } while(block<nBlocks);
  }

  void DMCPbyP::dmcWithBranching() {

    branchEngine->MaxAge=1;
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
        Mover->advanceWalkers(W.begin(),W.end());
        Eest = branchEngine->setWeights(W.begin(), W.end());
        ++step; ++CurrentStep;
        Estimators->accumulate(W);
        pop_acc += branchEngine->branch(CurrentStep,W);
        if(CurrentStep%100 == 0) updateWalkers();
      } while(step<nSteps);
      
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

  bool 
  DMCPbyP::put(xmlNodePtr q){
    return true;
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
