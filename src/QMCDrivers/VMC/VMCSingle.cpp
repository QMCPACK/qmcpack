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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/VMC/VMCSingle.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"

namespace qmcplusplus
{

/// Constructor.
VMCSingle::VMCSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool), Mover(0), UseDrift("yes")
{
  RootName = "vmc";
  QMCType ="VMPSingle";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  m_param.add(UseDrift,"useDrift","string");
  m_param.add(UseDrift,"usedrift","string");
  m_param.add(UseDrift,"use_drift","string");
  m_param.add(nTargetSamples,"targetWalkers","int");
  m_param.add(nTargetSamples,"targetwalkers","int");
  m_param.add(nTargetSamples,"target_walkers","int");
}

bool VMCSingle::run()
{
  resetRun();
  Mover->startRun(nBlocks,true);
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
  do
  {
    Mover->startBlock(nSteps);
    IndexType step = 0;
    do
    {
      ++step;
      ++CurrentStep;
      W.resetCollectables();
      Mover->advanceWalkers(W.begin(),W.end(),true); //step==nSteps);
      Estimators->accumulate(W);
      if(CurrentStep%updatePeriod==0)
        Mover->updateWalkers(W.begin(),W.end());
      if(CurrentStep%myPeriod4WalkerDump==0)
        W.saveEnsemble();
    }
    while(step<nSteps);
    Mover->stopBlock();
    nAcceptTot += Mover->nAccept;
    nRejectTot += Mover->nReject;
    ++block;
    recordBlock(block);
    ////periodically re-evaluate everything for pbyp
    //if(QMCDriverMode[QMC_UPDATE_MODE] && CurrentStep%100 == 0)
    //  Mover->updateWalkers(W.begin(),W.end());
  }
  while(block<nBlocks);
  Mover->stopRun();
  //finalize a qmc section
  return finalize(block);
}

void VMCSingle::resetRun()
{
  int samples_tot=W.getActiveWalkers()*nBlocks*nSteps*myComm->size();
  myPeriod4WalkerDump=(nTargetSamples>0)?samples_tot/nTargetSamples:Period4WalkerDump;
  if(myPeriod4WalkerDump==0 || QMCDriverMode[QMC_WARMUP])
    myPeriod4WalkerDump=(nBlocks+1)*nSteps;
  W.clearEnsemble();
  samples_tot=W.getActiveWalkers()*((nBlocks*nSteps)/myPeriod4WalkerDump);
  W.setNumSamples(samples_tot);
  if(Mover ==0)
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
    {
      app_log() << "  Update particle by particle " << endl;
      if(UseDrift == "yes")
        Mover=new VMCUpdatePbyPWithDrift(W,Psi,H,Random);
      else
        Mover=new VMCUpdatePbyP(W,Psi,H,Random);
    }
    else
    {
      app_log() << "  Update walker by walker " << endl;
      if(UseDrift == "yes")
        Mover=new VMCUpdateAllWithDrift(W,Psi,H,Random);
      else
        Mover=new VMCUpdateAll(W,Psi,H,Random);
    }
  }
  Mover->put(qmcNode);
  Mover->resetRun(branchEngine,Estimators);
  if(QMCDriverMode[QMC_UPDATE_MODE])
    Mover->initWalkersForPbyP(W.begin(),W.end());
  else
    Mover->initWalkers(W.begin(),W.end());
  app_log() << "  Samples are dumped at every " << myPeriod4WalkerDump << " step " << endl;
  app_log() << "  Total Sample Size =" << nTargetSamples
            << "\n  Sample size per node per thread = " << samples_tot << endl;
  app_log() << "  Warmup Steps " << nWarmupSteps << endl;
  //do a warmup run
  for(int prestep=0; prestep<nWarmupSteps; ++prestep)
    Mover->advanceWalkers(W.begin(),W.end(),true);
  nWarmupSteps=0;
}

bool
VMCSingle::put(xmlNodePtr q)
{
  //nothing to add
  return true;
}
}

/***************************************************************************
 * $RCSfile: VMCParticleByParticle.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCParticleByParticle.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/
