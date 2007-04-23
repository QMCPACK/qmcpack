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
#include "Message/OpenMP.h"
#define ENABLE_VMC_OMP_MASTER

namespace qmcplusplus { 

  /// Constructor.
  VMCSingleOMP::VMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
      HamiltonianPool& hpool):
    QMCDriver(w,psi,h),  CloneManager(hpool), UseDrift("yes") 
    { 
    RootName = "vmc";
    QMCType ="VMCSingleOMP";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
    m_param.add(UseDrift,"useDrift","string"); m_param.add(UseDrift,"usedrift","string");
  }
  
  /** execute the main body
   *
   * For eacth thread, cloned objects for branch, estimator and mover engines
   * are created at resetRun once during the life time of this object.
   * Master clause is used to collect estimators and perform file IO.
   */
#if defined(ENABLE_VMC_OMP_MASTER)
  bool VMCSingleOMP::run() 
  { 
    resetRun();

    //start the main estimator
    Estimators->start(nBlocks);

    //start a parallel region 
#pragma omp parallel 
    {
      int now=CurrentStep;
      int ip = omp_get_thread_num();
      bool pbyp=QMCDriverMode[QMC_UPDATE_MODE];
      //assign the iterators and resuse them
      MCWalkerConfiguration::iterator 
        wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
      //Creating active copy is not necessary
      //if(ip)
      //{
      //  if(wClones[ip]->getActiveWalkers() != (wPerNode[ip+1]-wPerNode[ip])) 
      //    wClones[ip]->createWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      //  wit=wClones[ip]->begin();
      //  wit_end=wClones[ip]->end();
      //}
      //else
      //{
      //  wit=W.begin();
      //  wit_end=W.begin()+wPerNode[1];
      //}
      if(pbyp)
        Movers[ip]->initWalkersForPbyP(wit,wit_end);
      else
        Movers[ip]->initWalkers(wit,wit_end);

      //disable recording function
      Movers[ip]->startRun(nBlocks,false);
      IndexType block = 0;
      do {
        Movers[ip]->startBlock(nSteps);
        IndexType step = 0;
        do 
        {
          ++step;
          ++now;
          Movers[ip]->advanceWalkers(wit,wit_end);
          Movers[ip]->accumulate(wit,wit_end);
        } while(step<nSteps);

        ++block;
        Movers[ip]->stopBlock();
#pragma omp barrier
#pragma omp master
        {
          Estimators->stopBlock(estimatorClones);
          recordBlock(block);
        }

        if(pbyp && now%100 == 0) 
        {
          Movers[ip]->updateWalkers(wit,wit_end);
        }
      } while(block<nBlocks);
    }

    CurrentStep+=nSteps*nBlocks;
    Estimators->stop(estimatorClones);
    //finalize a qmc section
    return finalize(nBlocks);
  }
#else /* !defined(ENABLE_VMC_OMP_MASTER) */
  bool VMCSingleOMP::run() 
  { 
    resetRun();

    //start the main estimator
    Estimators->start(nBlocks);

    for(int ip=0; ip<NumThreads; ip++) Movers[ip]->startRun(nBlocks,false);

    IndexType block = 0;
    do {
      //start a parallel region
#pragma omp parallel 
      {
        bool pbyp=QMCDriverMode[QMC_UPDATE_MODE];
        int ip = omp_get_thread_num();
        int now=CurrentStep;

        //assign the iterators and resuse them
        MCWalkerConfiguration::iterator 
          wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
        if(block ==0)
        {
          if(pbyp)
            Movers[ip]->initWalkersForPbyP(wit,wit_end);
          else
            Movers[ip]->initWalkers(wit,wit_end);
        }

        if(pbyp && now%100 == 99) Movers[ip]->updateWalkers(wit,wit_end);
        Movers[ip]->startBlock(nSteps);
        IndexType step = 0;
        do 
        {
          ++step;
          ++now;
          Movers[ip]->advanceWalkers(wit,wit_end);
          Movers[ip]->accumulate(wit,wit_end);
        } while(step<nSteps);
        Movers[ip]->stopBlock();
      }//end-of-parallel region

      CurrentStep+=nSteps;
      ++block;
      Estimators->stopBlock(estimatorClones);
      recordBlock(block);

    } while(block<nBlocks);

    Estimators->stop(estimatorClones);

    //finalize a qmc section
    return finalize(nBlocks);
  }
#endif /* ENABLE_VMC_OMP_MASTER */

  void VMCSingleOMP::resetRun() 
  {

    makeClones(W,Psi,H);

    if(Movers.empty()) 
    {
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
        if(ip) hClones[ip]->add2WalkerProperty(*wClones[ip]);
        estimatorClones[ip]= new ScalarEstimatorManager(*Estimators,*hClones[ip]);  
        estimatorClones[ip]->setCollectionMode(false);

        Rng[ip]=new RandomGenerator_t();
        Rng[ip]->init(OHMMS::Controller->mycontext()*NumThreads+ip,
            NumThreads*OHMMS::Controller->ncontexts(),-1);
        hClones[ip]->setRandomGenerator(Rng[ip]);

        branchClones[ip] = new BranchEngineType(*branchEngine);

        if(QMCDriverMode[QMC_UPDATE_MODE])
        {
          if(UseDrift == "yes")
            Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          else
            Movers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        }
        else
        {
          if(UseDrift == "yes")
            Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          else
            Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        }
      }
    }

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
