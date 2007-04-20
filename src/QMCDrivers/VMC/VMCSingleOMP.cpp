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
   * There exist one branchEngine and Estimators which is managed by the branchEngine
   * for a process. 
   */
  bool VMCSingleOMP::run() 
  { 
    resetRun();

    //start the main estimator
    Estimators->start(nBlocks);

    //start a parallel region
#pragma omp parallel  
    {
      int ip = omp_get_thread_num();
      if(QMCDriverMode[QMC_UPDATE_MODE])
        Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      else
        Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);

      //disable recording function
      Movers[ip]->startRun(nBlocks,false);
      IndexType block = 0;
      do {
        Movers[ip]->startBlock(nSteps);
        IndexType step = 0;
        do 
        {
          ++step;
          ++CurrentStep;
          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
          Movers[ip]->accumulate(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        } while(step<nSteps);

        ++block;
        Movers[ip]->stopBlock();
        //need a barrier???? 
//#pragma omp barrier
#pragma omp master
        {
          Estimators->stopBlock(estimatorClones);
          recordBlock(block);
        }

        if(QMCDriverMode[QMC_UPDATE_MODE] && CurrentStep%100 == 0) 
        {
          Movers[ip]->updateWalkers(W.begin()+wPerNode[ip], W.begin()+wPerNode[ip+1]);
        }
      } while(block<nBlocks);

      Movers[ip]->stopRun();
    }

    Estimators->stop(estimatorClones);

    //finalize a qmc section
    return finalize(nBlocks);
  }

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
        Rng[ip]=new RandomGenerator_t();
        Rng[ip]->init(ip,NumThreads,-1);
        hClones[ip]->setRandomGenerator(Rng[ip]);

        branchClones[ip] = new BranchEngineType(*branchEngine);

        if(QMCDriverMode[QMC_UPDATE_MODE])
        {
          if(UseDrift == "yes")
            Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          else
            Movers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
          //Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        }
        else
        {
          if(UseDrift == "yes")
            Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          else
            Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
          //Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
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
