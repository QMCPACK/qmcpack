//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/WFMC/WFMCSingleOMP.h"
#include "QMCDrivers/WFMC/WFMCUpdateAll.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
//#define ENABLE_VMC_OMP_MASTER

namespace qmcplusplus
{

/// Constructor.
WFMCSingleOMP::WFMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                             HamiltonianPool& hpool, WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool),  CloneManager(hpool),
  myWarmupSteps(0),UseDrift("yes"),reweight("yes")
{
  RootName = "wfmc";
  QMCType ="WFMCSingleOMP";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  m_param.add(UseDrift,"useDrift","string");
  m_param.add(UseDrift,"usedrift","string");
  m_param.add(myWarmupSteps,"warmupSteps","int");
  m_param.add(reweight,"reweight","string");
//     m_param.add(nTargetSamples,"targetWalkers","int");
}

bool WFMCSingleOMP::run()
{
  resetRun();
  //start the main estimator
  Estimators->start(nBlocks);
  ///Load a single walkers position into the walker.
  MCWalkerConfiguration Keeper(W);
  Keeper.createWalkers(W.begin(),W.end());
  MCWalkerConfiguration::iterator Kit=(Keeper.begin()), Kit_end(Keeper.end());
  int block=0;
  while ((Kit!=Kit_end)&&(nBlocks>block))
  {
    MCWalkerConfiguration::iterator Wit(W.begin()), Wit_end(W.end());
    while ((Wit!=Wit_end))
    {
      (*Wit)->R=(*Kit)->R;
      (*Wit)->G=(*Kit)->G;
      (*Wit)->L=(*Kit)->L;
      //(*Wit)->Drift=(*Kit)->Drift;
      (*Wit)->reset();
      (*Wit)->resetPropertyHistory();
      ++Wit;
    }
    #pragma omp parallel for
    for (int ip=0; ip<NumThreads; ++ip)
      Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    Wit=W.begin();
    while ((Wit!=Wit_end))
    {
      //       app_log()<<std::exp((*Wit)->Properties(LOGPSI))<< std::endl;
      (*Wit)->PropertyHistory[0][0]=std::exp((*Wit)->Properties(LOGPSI));
      ++Wit;
    }
    #pragma omp parallel
    {
      #pragma omp for
      for (int ip=0; ip<NumThreads; ++ip)
        Movers[ip]->startRun(nBlocks,false);
      #pragma omp for
      for (int ip=0; ip<NumThreads; ++ip)
      {
        //assign the iterators and resuse them
        MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
        Movers[ip]->startBlock(nSteps);
        for (int step=0; step<nSteps; ++step)
        {
          Movers[ip]->advanceWalkers(wit,wit_end,false);
          //      Movers[ip]->updateWalkers(wit,wit_end);
          //      wClones[ip]->saveEnsemble(wit,wit_end);
        }
        Movers[ip]->accumulate(wit,wit_end);
        Movers[ip]->stopBlock();
      }
      #pragma omp master
      {
        Estimators->stopBlock(estimatorClones);
        recordBlock(block+1);
        block++;
        ++Kit;
      }
    }
  }//end of parallel
  Estimators->stop(estimatorClones);
  //copy back the random states
  for (int ip=0; ip<NumThreads; ++ip)
    *(RandomNumberControl::Children[ip])=*(Rng[ip]);
  //finalize a qmc section
  return finalize(nBlocks);
}

void WFMCSingleOMP::resetRun()
{
  ///Set up a PropertyHistory for the energy to be recorded
  Eindex=0;
  MCWalkerConfiguration::iterator Cit(W.begin()), Cit_end(W.end());
  Eindex = (*Cit)->addPropertyHistory(nSteps);
  Cit++;
  while (Cit!=Cit_end)
  {
    (*Cit)->addPropertyHistory(nSteps);
    Cit++;
  }
  makeClones(W,Psi,H);
  //determine dump period for walkers
  int samples_tot=W.getActiveWalkers()*nBlocks*nSteps*myComm->size();
  myPeriod4WalkerDump=(nTargetSamples>0)?samples_tot/nTargetSamples:Period4WalkerDump;
  //fall back to the default
  if (myPeriod4WalkerDump==0)
    myPeriod4WalkerDump=Period4WalkerDump;
  if (QMCDriverMode[QMC_WARMUP])
    myPeriod4WalkerDump=nBlocks*nSteps;
  int samples_th=nTargetSamples/myComm->size()/NumThreads;
  for (int ip=0; ip<NumThreads; ++ip)
  {
    wClones[ip]->clearEnsemble();
    wClones[ip]->setNumSamples(samples_th);
  }
  app_log() << "  Samples are dumped at every " << myPeriod4WalkerDump << " step " << std::endl;
  app_log() << "  Total Sample Size =" << nTargetSamples
            << "\n  Sample size per node per thread = " << samples_th << std::endl;
  app_log() << "  Warmup Steps " << myWarmupSteps << std::endl;
  if (Movers.empty())
  {
    Movers.resize(NumThreads,0);
    branchClones.resize(NumThreads,0);
    estimatorClones.resize(NumThreads,0);
    Rng.resize(NumThreads,0);
    int nwtot=(W.getActiveWalkers()/NumThreads)*NumThreads;
    FairDivideLow(nwtot,NumThreads,wPerNode);
    app_log() << "  Initial partition of walkers ";
    copy(wPerNode.begin(),wPerNode.end(),std::ostream_iterator<int>(app_log()," "));
    app_log() << std::endl;
    #pragma omp parallel for
    for (int ip=0; ip<NumThreads; ++ip)
    {
      estimatorClones[ip]= new EstimatorManager(*Estimators);//,*hClones[ip]);
      estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
      estimatorClones[ip]->setCollectionMode(false);
      Rng[ip]=new RandomGenerator_t(*(RandomNumberControl::Children[ip]));
      hClones[ip]->setRandomGenerator(Rng[ip]);
      branchClones[ip] = new BranchEngineType(*branchEngine);
      if (reweight=="yes")
        Movers[ip]=new WFMCUpdateAllWithReweight(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip],nSteps,Eindex);
      else
        if (UseDrift == "yes")
          Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        else
          Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
    }
  }
  #pragma omp parallel  for
  for (int ip=0; ip<NumThreads; ++ip)
  {
    Movers[ip]->put(qmcNode);
    Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
    Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    for (int prestep=0; prestep<myWarmupSteps; ++prestep)
      Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
  }
  myWarmupSteps=0;
  //Used to debug and benchmark opnemp
  //#pragma omp parallel for
  //    for(int ip=0; ip<NumThreads; ip++)
  //    {
  //      Movers[ip]->benchMark(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],ip);
  //    }
}

bool
WFMCSingleOMP::put(xmlNodePtr q)
{
  //nothing to add
  return true;
}
}

