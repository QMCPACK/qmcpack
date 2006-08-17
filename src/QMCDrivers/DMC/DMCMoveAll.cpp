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
#include "QMCDrivers/DMC/DMCMoveAll.h"
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "QMCDrivers/DMC/DMCNonLocalUpdate.h"
#include "Message/Communicate.h"

namespace qmcplusplus {

  DMCMoveAll::DMCMoveAll(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h),Mover(0),
    BranchInterval(-1),KillNodeCrossing(0), NonLocalMoveIndex(-1),
    KillWalker("no"),NonLocalMove("no"), 
    Reconfiguration("no") { 
    RootName = "dmc";
    QMCType ="DMVMoveAll";
    m_param.add(KillWalker,"killnode","string");
    m_param.add(Reconfiguration,"reconfiguration","string");
    m_param.add(BranchInterval,"branchInterval","int"); 
    m_param.add(BranchInterval,"branch_interval","int");
    m_param.add(NonLocalMove,"nonlocalmoves","string");
  }

  DMCMoveAll::~DMCMoveAll() {
    if(Mover) delete Mover;
  }

  bool DMCMoveAll::put(xmlNodePtr cur){
    return true;
  }
  

  void DMCMoveAll::dmcWithBranching() {
    //m_oneover2tau = 0.5/Tau;
    //m_sqrttau = sqrt(Tau);
    
    RealType Eest = branchEngine->E_T;
    Mover->resetRun(branchEngine);
    Mover->MaxAge=3;

    bool checkNonLocalMove=(NonLocalMoveIndex>0);
    IndexType block = 0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;
    IndexType Population = W.getActiveWalkers();
    IndexType tPopulation = W.getActiveWalkers();
    do {
      IndexType step = 0;
      IndexType pop_acc=0; 

      Mover->startBlock();
      Estimators->startBlock();

      Mover->NonLocalMoveAccepted=0;
      do {
        pop_acc += W.getActiveWalkers();

        Mover->advanceWalkers(W.begin(),W.end());
        step++; CurrentStep++;
        Estimators->accumulate(W);
        Eest = branchEngine->CollectAndUpdate(W.getActiveWalkers(),Eest);
        branchEngine->branch(CurrentStep,W);
      } while(step<nSteps);


      nAccept = Mover->nAccept;
      nReject = Mover->nReject;
      RealType oneOverTotSteps=1.0/static_cast<RealType>(nAccept+nReject);
      if(checkNonLocalMove) {
        Estimators->setColumn(NonLocalMoveIndex, Mover->NonLocalMoveAccepted*oneOverTotSteps);
      }
      Estimators->stopBlock(static_cast<RealType>(nAccept)*oneOverTotSteps);

      nAcceptTot += nAccept;
      nRejectTot += nReject;
      
      Estimators->setColumn(PopIndex,static_cast<RealType>(pop_acc)/static_cast<RealType>(nSteps));
      Estimators->setColumn(EtrialIndex,Eest);

      Eest = Estimators->average(0);

      nAccept = 0; nReject = 0;
      block++;
      recordBlock(block);
    } while(block<nBlocks);
    
    //Need MPI-IO
    app_log() << "\t ratio = " << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot) << endl;
  }

  void DMCMoveAll::dmcWithReconfiguration() {
    Mover->MaxAge=0;
    IndexType block = 0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;
    RealType Eest = branchEngine->E_T;

    bool checkNonLocalMove=(NonLocalMoveIndex>0);
    //Mover->resetRun(branchEngine);
    do {
      IndexType step = 0;
      Mover->startBlock();
      Estimators->startBlock();
      do {
        int interval=0; 
        Mover->NonLocalMoveAccepted=0;
        do {
          Mover->advanceWalkers(W.begin(), W.end());
          ++interval;
          ++step; ++CurrentStep; 
        } while(interval<BranchInterval);

        if(checkNonLocalMove) {
          Estimators->setColumn(NonLocalMoveIndex, 
              static_cast<RealType>(Mover->NonLocalMoveAccepted)/
              static_cast<RealType>(W.getActiveWalkers()*BranchInterval));
        }

        Estimators->accumulate(W);
        int nwKept=branchEngine->branch(CurrentStep,W);
        Estimators->setColumn(PopIndex,nwKept);
        Eest = branchEngine->CollectAndUpdate(W.getActiveWalkers(),Eest);
      } while(step<nSteps);

      nAccept = Mover->nAccept;
      nReject = Mover->nReject;

      Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));

      nAcceptTot += nAccept;
      nRejectTot += nReject;
      
      Estimators->setColumn(EtrialIndex,Eest);

      Eest = Estimators->average(0);
      nAccept = 0; nReject = 0;
      block++;
      recordBlock(block);
      W.reset();
    } while(block<nBlocks);
    
    //Need MPI-IO
    app_log() << "\t ratio = " << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot) << endl;
  }
  /** Advance the walkers nblocks*nsteps timesteps. 
   *
   * For each block:
   * <ul>
   *  <li> Advance walkers for nsteps
   *  For each timestep:
   *   <ul>
   *   <li> Move all the particles of a walker.
   *   <li> Calculate the properties for the new walker configuration.
   *   <li> Accept/reject the new configuration.
   *   <li> Accumulate the estimators.
   *   <li> Update the trial energy \f$ E_T \f$
   *   <li> Branch the population of walkers (birth/death algorithm).
   *   </ul>
   * <li> Flush the estimators and print to file.
   * <li> Update the estimate of the local energy.
   * <li> (Optional) Print the ensemble of walker configurations.
   * </ul>
   * Default mode: Print the ensemble of walker configurations 
   * at the end of the run.
   */
  bool DMCMoveAll::run() { 

    bool fixW = (Reconfiguration == "yes");

    branchEngine->initWalkerController(Tau,fixW);
    KillNodeCrossing = (KillWalker == "yes");

    if(Mover ==0) {
      if(NonLocalMove == "yes") {
        app_log() << "  Non-local update is used." << endl;
        DMCNonLocalUpdate* nlocMover= new DMCNonLocalUpdate(W,Psi,H,Random);
        nlocMover->put(qmcNode);
        Mover=nlocMover;
        NonLocalMoveIndex=Estimators->addColumn("NonLocalMove");
      } else {
        if(KillNodeCrossing) {
          app_log() << "  Walkers will be killed if a node crossing is detected." << endl;
          Mover = new DMCUpdateAllWithKill(W,Psi,H,Random);
        } else {
          app_log() << "  Walkers will be kept even if a node crossing is detected." << endl;
          Mover = new DMCUpdateAllWithRejection(W,Psi,H,Random);
        }
      }
    }

    //add columns
    PopIndex = Estimators->addColumn("Population");
    EtrialIndex = Estimators->addColumn("Etrial");
    Estimators->reportHeader(AppendRun);

    Mover->resetRun(branchEngine);

    if(fixW)  {
      if(BranchInterval<0) {
        BranchInterval=nSteps;
        nSteps=1;
      }
      app_log() << "  DMC all-ptcl update with reconfigurations " << endl;
      app_log() << "    BranchInterval=" << BranchInterval << endl;
      app_log() << "    Steps         =" << nSteps << endl;
      app_log() << "    Blocks        =" << nBlocks << endl;
      dmcWithReconfiguration();
    } else {
      app_log() << "  DMC all-ptcl update with a fluctuating population" << endl;
      dmcWithBranching();
    }

    Estimators->finalize();
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
