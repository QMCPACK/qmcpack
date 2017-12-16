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
    
    
#include "QMCDrivers/DMC/DMCMoveAll.h"
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "QMCDrivers/DMC/DMCNonLocalUpdate.h"
#include "Estimators/DMCEnergyEstimator.h"

namespace qmcplusplus
{

DMCMoveAll::DMCMoveAll(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h),Mover(0),
  BranchInterval(-1),KillNodeCrossing(0), NonLocalMoveIndex(-1),
  KillWalker("no"),NonLocalMove("no"),
  Reconfiguration("no")
{
  RootName = "dmc";
  QMCType ="DMVMoveAll";
  m_param.add(KillWalker,"killnode","string");
  m_param.add(Reconfiguration,"reconfiguration","string");
  m_param.add(BranchInterval,"branchInterval","int");
  m_param.add(BranchInterval,"branch_interval","int");
  m_param.add(NonLocalMove,"nonlocalmove","string");
  m_param.add(NonLocalMove,"nonlocalmoves","string");
  //create a ScalarEstimator and add DMCEnergyEstimator
//     if (Estimators) delete Estimators;
  Estimators = new EstimatorManagerBase(H);
  Estimators->add(new DMCEnergyEstimator,"elocal");
}

DMCMoveAll::~DMCMoveAll()
{
  if(Mover)
    delete Mover;
}

bool DMCMoveAll::put(xmlNodePtr cur)
{
  return true;
}


bool DMCMoveAll::dmcWithBranching()
{
  //m_oneover2tau = 0.5/Tau;
  //m_sqrttau = sqrt(Tau);
  Mover->resetRun(branchEngine);
  Mover->MaxAge=3;
  bool checkNonLocalMove=(NonLocalMoveIndex>0);
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  IndexType Population = W.getActiveWalkers();
  IndexType tPopulation = W.getActiveWalkers();
  do
  {
    IndexType step = 0;
    IndexType pop_acc=0;
    Mover->startBlock();
    Estimators->startBlock();
    Mover->NonLocalMoveAccepted=0;
    do
    {
      pop_acc += W.getActiveWalkers();
      Mover->advanceWalkers(W.begin(),W.end());
      step++;
      CurrentStep++;
      Mover->setMultiplicity(W.begin(),W.end());
      Estimators->accumulate(W);
      branchEngine->branch(CurrentStep,W);
    }
    while(step<nSteps);
    nAccept = Mover->nAccept;
    nReject = Mover->nReject;
    RealType oneOverTotSteps=1.0/static_cast<RealType>(nAccept+nReject);
    if(checkNonLocalMove)
    {
      Estimators->setColumn(NonLocalMoveIndex, Mover->NonLocalMoveAccepted*oneOverTotSteps);
    }
    Estimators->stopBlock(static_cast<RealType>(nAccept)*oneOverTotSteps);
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    nAccept = 0;
    nReject = 0;
    block++;
    recordBlock(block);
  }
  while(block<nBlocks);
  //Need MPI-IO
  app_log() << "\t ratio = " << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot) << std::endl;
  return finalize(block);
}

bool DMCMoveAll::dmcWithReconfiguration()
{
  Mover->MaxAge=0;
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  bool checkNonLocalMove=(NonLocalMoveIndex>0);
  //Mover->resetRun(branchEngine);
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
        ++step;
        ++CurrentStep;
      }
      while(interval<BranchInterval);
      if(checkNonLocalMove)
      {
        Estimators->setColumn(NonLocalMoveIndex,
                              static_cast<RealType>(Mover->NonLocalMoveAccepted)/
                              static_cast<RealType>(W.getActiveWalkers()*BranchInterval));
      }
      Estimators->accumulate(W);
      branchEngine->branch(CurrentStep,W);
    }
    while(step<nSteps);
    nAccept = Mover->nAccept;
    nReject = Mover->nReject;
    Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    nAccept = 0;
    nReject = 0;
    block++;
    recordBlock(block);
    W.reset();
  }
  while(block<nBlocks);
  //Need MPI-IO
  app_log() << "\t ratio = " << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot) << std::endl;
  return finalize(block);
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
bool DMCMoveAll::run()
{
  bool fixW = (Reconfiguration == "yes");
  branchEngine->initWalkerController(Tau,fixW);
  KillNodeCrossing = (KillWalker == "yes");
  if(Mover ==0)
  {
    if(NonLocalMove == "yes")
    {
      app_log() << "  Non-local update is used." << std::endl;
      DMCNonLocalUpdate* nlocMover= new DMCNonLocalUpdate(W,Psi,H,Random);
      nlocMover->put(qmcNode);
      Mover=nlocMover;
      NonLocalMoveIndex=Estimators->addColumn("NonLocalMove");
    }
    else
    {
      if(KillNodeCrossing)
      {
        app_log() << "  Walkers will be killed if a node crossing is detected." << std::endl;
        Mover = new DMCUpdateAllWithKill(W,Psi,H,Random);
      }
      else
      {
        app_log() << "  Walkers will be kept even if a node crossing is detected." << std::endl;
        Mover = new DMCUpdateAllWithRejection(W,Psi,H,Random);
      }
    }
  }
  Estimators->reportHeader(AppendRun);
  Mover->resetRun(branchEngine);
  bool success=false;
  if(fixW)
  {
    if(BranchInterval<0)
    {
      BranchInterval=nSteps;
      nSteps=1;
    }
    app_log() << "  DMC all-ptcl update with reconfigurations " << std::endl;
    app_log() << "    BranchInterval=" << BranchInterval << std::endl;
    app_log() << "    Steps         =" << nSteps << std::endl;
    app_log() << "    Blocks        =" << nBlocks << std::endl;
    success = dmcWithReconfiguration();
  }
  else
  {
    app_log() << "  DMC all-ptcl update with a fluctuating population" << std::endl;
    success = dmcWithBranching();
  }
  return success;
}
}
