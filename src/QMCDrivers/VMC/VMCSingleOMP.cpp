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
  
  void VMCSingleOMP::resetRun() 
  {

    makeClones(W,Psi,H);

    if(Movers.empty()) 
    {
      Movers.resize(NumThreads,0);
      branchClones.resize(NumThreads,0);
      Rng.resize(NumThreads,0);
      FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);

      app_log() << "  Initial partition of walkers ";
      std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(app_log()," "));
      app_log() << endl;

#pragma omp parallel  
      {
        int ip = omp_get_thread_num();
        if(ip) {
          hClones[ip]->add2WalkerProperty(*wClones[ip]);
        }
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
          Movers[ip]->resetRun(branchClones[ip]);
          Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        }
        else
        {
          if(UseDrift == "yes")
            Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          else
            Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
          Movers[ip]->resetRun(branchClones[ip]);
          Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        }
      }
    }

    for(int ip=0;ip<NumThreads; ip++) Movers[ip]->put(qmcNode);

    //Used to debug and benchmark opnemp
//#pragma omp parallel for
//    for(int ip=0; ip<NumThreads; ip++)
//    {
//      Movers[ip]->benchMark(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],ip);
//    }
  }

  bool VMCSingleOMP::run() 
  { 

    resetRun();

    Estimators->reportHeader(AppendRun);
    Estimators->reset();

    IndexType block = 0;
    do {
      Estimators->startBlock();
      for(int ip=0; ip<NumThreads; ip++) Movers[ip]->startBlock();

      IndexType step = 0;
      do 
      {
#pragma omp parallel
        {
          int ip = omp_get_thread_num();
          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        }
        ++step;
        ++CurrentStep;
        Estimators->accumulate(W);
      } while(step<nSteps);

      IndexType nAcceptTot=0;
      IndexType nRejectTot=0;
      for(int ip=0; ip<NumThreads; ip++) {
        nAcceptTot+=Movers[ip]->nAccept;
        nRejectTot+=Movers[ip]->nReject;
      }

      Estimators->stopBlock(static_cast<RealType>(nAcceptTot)/static_cast<RealType>(nAcceptTot+nRejectTot));
      ++block;

      //record the current configuration
      recordBlock(block);

      if(QMCDriverMode[QMC_UPDATE_MODE] && CurrentStep%100 == 0) 
      {
#pragma omp parallel  
        {
          int ip = omp_get_thread_num();
          Movers[ip]->updateWalkers(W.begin()+wPerNode[ip], W.begin()+wPerNode[ip+1]);
        }
      }
    } while(block<nBlocks);

    //finalize a qmc section
    return finalize(block);
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
