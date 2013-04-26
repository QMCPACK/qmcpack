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
#include "QMCDrivers/DMC/DMCOMP.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
//#include "QMCDrivers/DMC/DMCNonLocalUpdate.h"
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Utilities/Timer.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{

/// Constructor.
DMCOMP::DMCOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool)
  : QMCDriver(w,psi,h,ppool), CloneManager(hpool)
  , KillNodeCrossing(0) ,Reconfiguration("no"), BenchMarkRun("no"), UseFastGrad("yes")
  , BranchInterval(-1),mover_MaxAge(-1)
{
  RootName = "dmc";
  QMCType ="DMCOMP";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  m_param.add(KillWalker,"killnode","string");
  m_param.add(BenchMarkRun,"benchmark","string");
  m_param.add(Reconfiguration,"reconfiguration","string");
  m_param.add(BranchInterval,"branchInterval","string");
  m_param.add(NonLocalMove,"nonlocalmove","string");
  m_param.add(NonLocalMove,"nonlocalmoves","string");
  m_param.add(mover_MaxAge,"MaxAge","double");
  m_param.add(UseFastGrad,"fastgrad", "string");
  //DMC overwrites ConstPopulation
  ConstPopulation=false;
}

void DMCOMP::resetComponents(xmlNodePtr cur)
{
  qmcNode=cur;
  m_param.put(cur);
  put(cur);
  //app_log()<<"DMCOMP::resetComponents"<<endl;
  Estimators->reset();
  branchEngine->resetRun(cur);
  branchEngine->checkParameters(W);
  //delete Movers[0];
  for(int ip=0; ip<NumThreads; ++ip)
  {
    delete Movers[ip];
    delete estimatorClones[ip];
    delete branchClones[ip];
    estimatorClones[ip]= new EstimatorManager(*Estimators);
    estimatorClones[ip]->setCollectionMode(false);
    branchClones[ip] = new BranchEngineType(*branchEngine);
  }
#if !defined(BGP_BUG)
  #pragma omp parallel for
#endif
  for(int ip=0; ip<NumThreads; ++ip)
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
    {
      if(UseFastGrad == "yes")
        Movers[ip] = new DMCUpdatePbyPWithRejectionFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
      else
        Movers[ip] = new DMCUpdatePbyPWithRejection(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
      Movers[ip]->put(cur);
      Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    }
    else
    {
      if(KillNodeCrossing)
        Movers[ip] = new DMCUpdateAllWithKill(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
      else
        Movers[ip] = new DMCUpdateAllWithRejection(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
      Movers[ip]->put(cur);
      Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    }
  }
}


void DMCOMP::resetUpdateEngines()
{
  ReportEngine PRE("DMCOMP","resetUpdateEngines");
  bool fixW = (Reconfiguration == "yes");
  makeClones(W,Psi,H);
  Timer init_timer;
  if(Movers.empty())
  {
    W.loadEnsemble(wClones);
    branchEngine->initWalkerController(W,fixW,false);
    //if(QMCDriverMode[QMC_UPDATE_MODE]) W.clearAuxDataSet();
    Movers.resize(NumThreads,0);
    branchClones.resize(NumThreads,0);
    Rng.resize(NumThreads,0);
    estimatorClones.resize(NumThreads,0);
    FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
    {
      //log file
      ostringstream o;
      o << "  Initial partition of walkers on a node: ";
      std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(o," "));
      o << "\n";
      if(QMCDriverMode[QMC_UPDATE_MODE])
        o << "  Updates by particle-by-particle moves";
      else
        o << "  Updates by walker moves";
      if(UseFastGrad == "yes")
        o << " using fast gradient version ";
      else
        o << " using full-ratio version ";
      if(KillNodeCrossing)
        o << "\n  Walkers are killed when a node crossing is detected";
      else
        o << "\n  DMC moves are rejected when a node crossing is detected";
      app_log() << o.str() << endl;
    }
#if !defined(BGP_BUG)
    #pragma omp parallel for
#endif
    for(int ip=0; ip<NumThreads; ++ip)
    {
      estimatorClones[ip]= new EstimatorManager(*Estimators);
      estimatorClones[ip]->setCollectionMode(false);
      Rng[ip]=new RandomGenerator_t(*RandomNumberControl::Children[ip]);
      hClones[ip]->setRandomGenerator(Rng[ip]);
      branchClones[ip] = new BranchEngineType(*branchEngine);
      if(QMCDriverMode[QMC_UPDATE_MODE])
      {
        if(UseFastGrad == "yes")
          Movers[ip] = new DMCUpdatePbyPWithRejectionFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        else
          Movers[ip] = new DMCUpdatePbyPWithRejection(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        Movers[ip]->put(qmcNode);
        Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      }
      else
      {
        if(KillNodeCrossing)
          Movers[ip] = new DMCUpdateAllWithKill(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        else
          Movers[ip] = new DMCUpdateAllWithRejection(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        Movers[ip]->put(qmcNode);
        Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      }
    }
  }
  branchEngine->checkParameters(W);
  int mxage=mover_MaxAge;
  if(fixW)
  {
    if(BranchInterval<0)
      BranchInterval=nSteps;
    mxage=(mover_MaxAge<0)?0:mover_MaxAge;
    for(int ip=0; ip<Movers.size(); ++ip)
      Movers[ip]->MaxAge=mxage;
  }
  else
  {
    if(BranchInterval<0)
      BranchInterval=1;
    int miage=(QMCDriverMode[QMC_UPDATE_MODE])?1:5;
    mxage=(mover_MaxAge<0)?miage:mover_MaxAge;
    for(int i=0; i<NumThreads; ++i)
      Movers[i]->MaxAge=mxage;
  }
  {
    ostringstream o;
    if(fixW)
      o << "  Fixed population using reconfiguration method\n";
    else
      o << "  Fluctuating population\n";
    o << "  Persisent walkers are killed after " << mxage << " MC sweeps\n";
    o << "  BranchInterval = " << BranchInterval << "\n";
    o << "  Steps per block = " << nSteps << "\n";
    o << "  Number of blocks = " << nBlocks << "\n";
    app_log() << o.str() << endl;
  }
  app_log() << "  DMC Engine Initialization = " << init_timer.elapsed() << " secs " << endl;
}

bool DMCOMP::run()
{
  bool variablePop = (Reconfiguration == "no");
  resetUpdateEngines();
  //estimator does not need to collect data
  Estimators->setCollectionMode(true);
  Estimators->start(nBlocks);
  for(int ip=0; ip<NumThreads; ip++)
    Movers[ip]->startRun(nBlocks,false);
  Timer myclock;
  IndexType block = 0;
  IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
  do // block
  {
    Estimators->startBlock(nSteps);
    for(int ip=0; ip<NumThreads; ip++)
      Movers[ip]->startBlock(nSteps);
    IndexType step = 0;
    for(IndexType step=0; step< nSteps; ++step, CurrentStep+=BranchInterval)
    {
//         if(storeConfigs && (CurrentStep%storeConfigs == 0)) {
//           ForwardWalkingHistory.storeConfigsForForwardWalking(W);
//           W.resetWalkerParents();
//         }
      #pragma omp parallel
      {
        int ip=omp_get_thread_num();
        int now=CurrentStep;
        MCWalkerConfiguration::iterator
        wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
        for(int interval = 0; interval<BranchInterval-1; ++interval,++now)
          Movers[ip]->advanceWalkers(wit,wit_end,false);
        wClones[ip]->resetCollectables();
        Movers[ip]->advanceWalkers(wit,wit_end,false);
        Movers[ip]->setMultiplicity(wit,wit_end);
        if(QMCDriverMode[QMC_UPDATE_MODE] && now%updatePeriod == 0)
          Movers[ip]->updateWalkers(wit, wit_end);
      }//#pragma omp parallel
      //Collectables are weighted but not yet normalized
      if(W.Collectables.size())
      {
        // only when collectable is not empty, need to generalize for W != wClones[0]
        for(int ip=1; ip<NumThreads; ++ip)
          W.Collectables += wClones[ip]->Collectables;
      }
      branchEngine->branch(CurrentStep, W, branchClones);
//         if(storeConfigs && (CurrentStep%storeConfigs == 0)) {
//           ForwardWalkingHistory.storeConfigsForForwardWalking(W);
//           W.resetWalkerParents();
//         }
      if(variablePop)
        FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
    }
//       branchEngine->debugFWconfig();
    Estimators->stopBlock(acceptRatio());
    block++;
    if(DumpConfig &&block%Period4CheckPoint == 0)
    {
      for(int ip=0; ip<NumThreads; ip++)
        *(RandomNumberControl::Children[ip])=*(Rng[ip]);
    }
    recordBlock(block);
  }
  while(block<nBlocks && myclock.elapsed()<MaxCPUSecs);
  //for(int ip=0; ip<NumThreads; ip++) Movers[ip]->stopRun();
  for(int ip=0; ip<NumThreads; ip++)
    *(RandomNumberControl::Children[ip])=*(Rng[ip]);
  Estimators->stop();
  return finalize(block);
}

void DMCOMP::benchMark()
{
  //set the collection mode for the estimator
  Estimators->setCollectionMode(true);
  IndexType PopIndex = Estimators->addProperty("Population");
  IndexType EtrialIndex = Estimators->addProperty("Etrial");
  //Estimators->reportHeader(AppendRun);
  //Estimators->reset();
  IndexType block = 0;
  RealType Eest = branchEngine->getEref();
  //resetRun();
  for(int ip=0; ip<NumThreads; ip++)
  {
    char fname[16];
    sprintf(fname,"test.%i",ip);
    ofstream fout(fname);
  }
  for(int istep=0; istep<nSteps; istep++)
  {
    FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
    #pragma omp parallel
    {
      int ip = omp_get_thread_num();
      Movers[ip]->benchMark(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],ip);
    }
  }
}

bool
DMCOMP::put(xmlNodePtr q)
{
  //nothing to do
  return true;
}
}

/***************************************************************************
 * $RCSfile: DMCOMP.cpp,v $   $Author: jnkim $
 * $Revision: 1620 $   $Date: 2007-01-14 18:12:23 -0600 (Sun, 14 Jan 2007) $
 * $Id: DMCOMP.cpp 1620 2007-01-15 00:12:23Z jnkim $
 ***************************************************************************/
