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
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
//#define ENABLE_VMC_OMP_MASTER

namespace qmcplusplus { 

  /// Constructor.
  VMCSingleOMP::VMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
      HamiltonianPool& hpool):
    QMCDriver(w,psi,h),  CloneManager(hpool), 
    myWarmupSteps(0),UseDrift("yes") ,weightLength(0), reweight("no"), Eindex(-1)
  { 
    RootName = "vmc";
    QMCType ="VMCSingleOMP";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
    QMCDriverMode.set(QMC_WARMUP,0);
    m_param.add(UseDrift,"useDrift","string"); m_param.add(UseDrift,"usedrift","string");
    m_param.add(reweight,"reweight","string");
    m_param.add(weightLength,"weightlength","int");
    m_param.add(myWarmupSteps,"warmupSteps","int");
    m_param.add(nTargetSamples,"targetWalkers","int");
  }

  bool VMCSingleOMP::run() 
  { 
    resetRun();

    //start the main estimator
    Estimators->start(nBlocks);

    for(int ip=0; ip<NumThreads; ++ip) Movers[ip]->startRun(nBlocks,false);

    for(int block=0;block<nBlocks; ++block)
    {
#pragma omp parallel for 
      for(int ip=0; ip<NumThreads; ++ip)
      {
        //IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
        IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:0;
        //assign the iterators and resuse them
        MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);

        Movers[ip]->startBlock(nSteps);
        int now_loc=CurrentStep;
        for(int step=0; step<nSteps;++step)
        {
          wClones[ip]->resetCollectables();
          Movers[ip]->advanceWalkers(wit,wit_end,false);
          Movers[ip]->accumulate(wit,wit_end);
          ++now_loc;
          if(updatePeriod&& now_loc%updatePeriod==0) Movers[ip]->updateWalkers(wit,wit_end);
          if(Period4WalkerDump&& now_loc%myPeriod4WalkerDump==0) wClones[ip]->saveEnsemble(wit,wit_end);
        } 
        Movers[ip]->stopBlock(false);
      }//end-of-parallel for

      CurrentStep+=nSteps;
      Estimators->stopBlock(estimatorClones);
      //recordBlock(block+1);
    }//block

    Estimators->stop(estimatorClones);

    //copy back the random states
    for(int ip=0; ip<NumThreads; ++ip) 
      *(RandomNumberControl::Children[ip])=*(Rng[ip]);

    //finalize a qmc section
    return finalize(nBlocks);
  }

  void VMCSingleOMP::resetRun() 
  {
    
    ///Set up a PropertyHistory for the energy to be recorded
    if ((reweight=="yes") & (Eindex<0))
    {
      MCWalkerConfiguration::iterator Cit(W.begin()), Cit_end(W.end());
      Eindex = (*Cit)->addPropertyHistory(weightLength);
      while(Cit!=Cit_end){
	(*Cit)->addPropertyHistory(weightLength);
	Cit++;
      }
    }
    makeClones(W,Psi,H);

    //determine dump period for walkers
    int samples_tot=W.getActiveWalkers()*nBlocks*nSteps*myComm->size();
    myPeriod4WalkerDump=(nTargetSamples>0)?samples_tot/nTargetSamples:Period4WalkerDump;
    //fall back to the default
    if(myPeriod4WalkerDump==0) myPeriod4WalkerDump=Period4WalkerDump;
    if(QMCDriverMode[QMC_WARMUP]) myPeriod4WalkerDump=nBlocks*nSteps;
    int samples_th=nTargetSamples/myComm->size()/NumThreads;
    app_log() << "  Samples are dumped at every " << myPeriod4WalkerDump << " step " << endl;
    app_log() << "  Total Sample Size =" << nTargetSamples
      << "\n  Sample size per node per thread = " << samples_th << endl;
    app_log() << "  Warmup Steps " << myWarmupSteps << endl;

    if(Movers.empty()) 
    {
      Movers.resize(NumThreads,0);
      branchClones.resize(NumThreads,0);
      estimatorClones.resize(NumThreads,0);
      Rng.resize(NumThreads,0);
      int nwtot=(W.getActiveWalkers()/NumThreads)*NumThreads;
      FairDivideLow(nwtot,NumThreads,wPerNode);

      app_log() << "  Initial partition of walkers ";
      std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(app_log()," "));
      app_log() << endl;

#pragma omp parallel for
      for(int ip=0; ip<NumThreads; ++ip)
      {
        estimatorClones[ip]= new EstimatorManager(*Estimators);//,*hClones[ip]);  
        estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
        estimatorClones[ip]->setCollectionMode(false);
        Rng[ip]=new RandomGenerator_t(*(RandomNumberControl::Children[ip]));
        hClones[ip]->setRandomGenerator(Rng[ip]);

        branchClones[ip] = new BranchEngineType(*branchEngine);

        if(reweight=="yes")
        {
          app_log() << "  WFMCUpdateAllWithReweight"<<endl;
          Movers[ip]=new WFMCUpdateAllWithReweight(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip],weightLength,Eindex);
        }
        else if (reweight=="psi")
        {
          app_log() << "  Sampling Psi to increase number of walkers near nodes"<<endl;
          if (QMCDriverMode[QMC_UPDATE_MODE]) Movers[ip]=new VMCUpdatePbyPSamplePsi(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
          else Movers[ip]=new VMCUpdateAllSamplePsi(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        else if(QMCDriverMode[QMC_UPDATE_MODE])
        {
          app_log() << "  PbyP moves"<<endl;
          if(UseDrift == "yes")
            Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          else
            Movers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        }
        else
        {
          app_log() << "  Walker moves"<<endl;
          if(UseDrift == "yes")
            Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          else
            Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        }
      }
    }

#pragma omp parallel  for
    for(int ip=0; ip<NumThreads; ++ip)
    {
      Movers[ip]->put(qmcNode);
      Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);

      if(QMCDriverMode[QMC_UPDATE_MODE])
        Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      else
        Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);

      for(int prestep=0; prestep<myWarmupSteps; ++prestep)
        Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true); 

      if(myWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
        Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]); 

#pragma omp critical
      {
        wClones[ip]->clearEnsemble();
        wClones[ip]->setNumSamples(samples_th);
      }
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
  VMCSingleOMP::put(xmlNodePtr q){
    //nothing to add
    return true;
  }
}

/***************************************************************************
 * $RCSfile: VMCSingleOMP.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCSingleOMP.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $ 
 ***************************************************************************/
