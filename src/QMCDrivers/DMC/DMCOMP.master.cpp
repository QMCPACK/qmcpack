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
#include "QMCDrivers/DMC/DMCNonLocalUpdate.h"
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

/// Constructor.
DMCOMP::DMCOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
               HamiltonianPool& hpool):
  QMCDriver(w,psi,h), CloneManager(hpool),
  KillNodeCrossing(0),
  Reconfiguration("no"), BenchMarkRun("no"), BranchInterval(-1)
{
  RootName = "dummy";
  QMCType ="DMCOMP";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_MULTIPLE,1);
  m_param.add(KillWalker,"killnode","string");
  m_param.add(BenchMarkRun,"benchmark","string");
  m_param.add(Reconfiguration,"reconfiguration","string");
  m_param.add(BranchInterval,"branchInterval","string");
}

void DMCOMP::resetUpdateEngines()
{
  bool fixW = (Reconfiguration == "yes");
  makeClones(W,Psi,H);
  if(Movers.empty())
  {
    branchEngine->initWalkerController(Tau,fixW);
    //if(QMCDriverMode[QMC_UPDATE_MODE]) W.clearAuxDataSet();
    Movers.resize(NumThreads,0);
    branchClones.resize(NumThreads,0);
    estimatorClones.resize(NumThreads,0);
    Rng.resize(NumThreads,0);
    FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
    app_log() << "  Initial partition of walkers ";
    std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(app_log()," "));
    app_log() << endl;
    #pragma omp parallel
    {
      int ip = omp_get_thread_num();
      if(ip)
        hClones[ip]->add2WalkerProperty(*wClones[ip]);
      estimatorClones[ip]= new ScalarEstimatorManager(*Estimators,*hClones[ip]);
      Rng[ip]=new RandomGenerator_t();
      Rng[ip]->init(ip,NumThreads,-1);
      hClones[ip]->setRandomGenerator(Rng[ip]);
      branchClones[ip] = new BranchEngineType(*branchEngine);
      if(QMCDriverMode[QMC_UPDATE_MODE])
      {
        if(NonLocalMove == "yes")
        {
          DMCNonLocalUpdatePbyP* nlocMover=
            new DMCNonLocalUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
          nlocMover->put(qmcNode);
          Movers[ip]=nlocMover;
        }
        else
        {
          Movers[ip]= new DMCUpdatePbyPWithRejection(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        //Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      }
      else
      {
        if(NonLocalMove == "yes")
        {
          app_log() << "  Non-local update is used." << endl;
          DMCNonLocalUpdate* nlocMover= new DMCNonLocalUpdate(W,Psi,H,Random);
          nlocMover->put(qmcNode);
          Movers[ip]=nlocMover;
        }
        else
        {
          if(KillNodeCrossing)
          {
            Movers[ip] = new DMCUpdateAllWithKill(W,Psi,H,Random);
          }
          else
          {
            Movers[ip] = new DMCUpdateAllWithRejection(W,Psi,H,Random);
          }
        }
        Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        //Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      }
    }
  }
  if(fixW)
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
      app_log() << "  DMCOMP PbyP Update with reconfigurations" << endl;
    else
      app_log() << "  DMCOMP walker Update with reconfigurations" << endl;
    for(int ip=0; ip<Movers.size(); ip++)
      Movers[ip]->MaxAge=0;
    if(BranchInterval<0)
    {
      BranchInterval=nSteps;
    }
  }
  else
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
    {
      app_log() << "  DMCOMP PbyP Update with a fluctuating population" << endl;
      for(int ip=0; ip<Movers.size(); ip++)
        Movers[ip]->MaxAge=1;
    }
    else
    {
      app_log() << "  DMCOMP walker Update with a fluctuating population" << endl;
      for(int ip=0; ip<Movers.size(); ip++)
        Movers[ip]->MaxAge=3;
    }
    if(BranchInterval<0)
      BranchInterval=1;
  }
  app_log() << "  BranchInterval = " << BranchInterval << endl;
  app_log() << "  Steps per block = " << nSteps << endl;
  app_log() << "  Number of blocks = " << nBlocks << endl;
}

bool DMCOMP::run()
{
  bool variablePop = (Reconfiguration == "no");
  resetUpdateEngines();
  Estimators->start(nBlocks);
  #pragma omp parallel
  {
    int ip = omp_get_thread_num();
    int np = omp_get_num_threads();
    vector<int> wpart(np+1);
    FairDivideLow(W.getActiveWalkers(),np,wpart);
    int firstW=wpart[ip];
    int lastW=wpart[ip+1];
    //char fname[16];
    //sprintf(fname,"debug.%d",ip);
    //ofstream fout(fname);
    if(QMCDriverMode[QMC_UPDATE_MODE])
      Movers[ip]->initWalkersForPbyP(W.begin()+firstW,W.begin()+lastW);
    else
      Movers[ip]->initWalkers(W.begin()+firstW,W.begin()+lastW);
    Movers[ip]->startRun(nBlocks,false);
    IndexType block = 0;
    do // block
    {
      Movers[ip]->startBlock(nSteps);
      IndexType step = 0;
      do  //step
      {
        IndexType interval = 0;
        IndexType bi=BranchInterval;
        do // interval
        {
          ++interval;
          Movers[ip]->advanceWalkers(W.begin()+firstW,W.begin()+lastW);
        }
        while(interval<bi);
        Movers[ip]->setMultiplicity(W.begin()+firstW,W.begin()+lastW);
        #pragma omp barrier
        #pragma omp master
        {
          //branchEngine->branch(CurrentStep,W, branchClones);
          branchEngine->branch(CurrentStep,W);
          if(storeConfigs && (CurrentStep%storeConfigs == 0))
          {
            branchEngine->storeConfigsForForwardWalking(W);
            W.resetWalkerParents();
          }
        }
        #pragma omp barrier
        branchClones[ip]->E_T=branchEngine->E_T;
        if(variablePop)
        {
          FairDivideLow(W.getActiveWalkers(),np,wpart);
          firstW=wpart[ip];
          lastW=wpart[ip+1];
        }
        ++step;
//fout <<  step << " " << firstW << " " << lastW << " " << W.getActiveWalkers() << " " << branchClones[ip]->E_T << endl;
        CurrentStep+=BranchInterval;
      }
      while(step<nSteps);
      block++;
      #pragma omp master
      {
        Estimators->stopBlock(Movers[0]->acceptRatio());
        recordBlock(block);
      }
      if(QMCDriverMode[QMC_UPDATE_MODE] && CurrentStep%100 == 0)
      {
        Movers[ip]->updateWalkers(W.begin()+firstW, W.begin()+lastW);
      }
    }
    while(block<nBlocks);
    Movers[ip]->stopRun();
  }
  Estimators->stop();
  return finalize(nBlocks);
}

void DMCOMP::benchMark()
{
  //set the collection mode for the estimator
  Estimators->setCollectionMode(branchEngine->SwapMode);
  IndexType PopIndex = Estimators->addColumn("Population");
  IndexType EtrialIndex = Estimators->addColumn("Etrial");
  //Estimators->reportHeader(AppendRun);
  //Estimators->reset();
  IndexType block = 0;
  RealType Eest = branchEngine->E_T;
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

//  void DMCOMP::dmcWithBranching() {
//
//    RealType Eest = branchEngine->E_T;
//
//    resetRun();
//
//
//    IndexType block = 0;
//    do {//start a block
//      for(int ip=0; ip<NumThreads; ip++) {
//        Movers[ip]->startBlock();
//      }
//
//      FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
//      Estimators->startBlock();
//      IndexType step = 0;
//      IndexType pop_acc=0;
//      do {
//#pragma omp parallel
//        {
//          int ip = omp_get_thread_num();
//          Movers[ip]->resetEtrial(Eest);
//          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip], W.begin()+wPerNode[ip+1]);
//        }
//
//        ++step; ++CurrentStep;
//        Estimators->accumulate(W);
//
//        //int cur_pop = branchEngine->branch(CurrentStep,W, branchClones);
//        //pop_acc += cur_pop;
//        //Eest = branchEngine->CollectAndUpdate(cur_pop, Eest);
//
//        FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
//
//        for(int ip=0; ip<NumThreads; ip++) Movers[ip]->resetEtrial(Eest);
//
//        if(CurrentStep%100 == 0) {
//#pragma omp parallel
//          {
//            int ip = omp_get_thread_num();
//            Movers[ip]->updateWalkers(W.begin()+wPerNode[ip], W.begin()+wPerNode[ip+1]);
//          }
//        }
//      } while(step<nSteps);
//
//      Estimators->stopBlock(acceptRatio());
//
//      //update estimator
//      //Estimators->setColumn(PopIndex,static_cast<RealType>(pop_acc)/static_cast<RealType>(nSteps));
//      //Estimators->setColumn(EtrialIndex,Eest);
//      //Eest = Estimators->average(0);
//      //RealType totmoves=1.0/static_cast<RealType>(step*W.getActiveWalkers());
//
//      //Need MPI-IO
//      //app_log()
//      //  << setw(4) << block
//      //  << setw(20) << static_cast<RealType>(nAllRejected)*totmoves
//      //  << setw(20) << static_cast<RealType>(nNodeCrossing)*totmoves << endl;
//      block++;
//
//      recordBlock(block);
//
//    } while(block<nBlocks);
//  }
//


}

/***************************************************************************
 * $RCSfile: DMCOMP.cpp,v $   $Author: jnkim $
 * $Revision: 1620 $   $Date: 2007-01-14 18:12:23 -0600 (Sun, 14 Jan 2007) $
 * $Id: DMCOMP.cpp 1620 2007-01-15 00:12:23Z jnkim $
 ***************************************************************************/
