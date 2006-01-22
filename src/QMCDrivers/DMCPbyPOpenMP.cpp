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
#include "QMCDrivers/DMCPbyPOpenMP.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"

namespace qmcplusplus { 

  /// Constructor.
  DMCPbyPOpenMP::DMCPbyPOpenMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h){ 
    RootName = "dummy";
    QMCType ="dummy";
    NumThreads=omp_get_max_threads();
    wPerNode.resize(NumThreads+1,0);
  }

  void DMCPbyPOpenMP::makeClones(HamiltonianPool& hpool, int np) {

    //np is ignored but will see if ever needed
    wPerNode.resize(NumThreads+1,0);

    if(wClones.size()) {
      ERRORMSG("Cannot make clones again. Use clones")
      return;
    }

    wClones.resize(NumThreads,0);
    psiClones.resize(NumThreads,0);
    hClones.resize(NumThreads,0);
    Rng.resize(NumThreads,0);

    wClones[0]=&W;
    psiClones[0]=&Psi;
    hClones[0]=&H;

    if(NumThreads == 1) {
      WARNMSG("Using a single thread with DMCPbyPOpenMP.")
      return;
    }

    hpool.clone(W,Psi,H,wClones,psiClones,hClones);

  }

  void DMCPbyPOpenMP::resetRun() {

    if(Movers.empty()) {
      Movers.resize(NumThreads,0);
#pragma omp parallel  
      {
        int ip = omp_get_thread_num();
        Rng[ip]=new RandomGenerator_t();
        Rng[ip]->init(ip,NumThreads,-1);
        Movers[ip] = new DMCPbyPUpdate(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
        Movers[ip]->resetRun(new BranchEngineType(*branchEngine));
      }
    }
  }
  
  bool DMCPbyPOpenMP::run() { 
    
    //set the collection mode for the estimator
    Estimators->setCollectionMode(branchEngine->SwapMode);

    IndexType PopIndex = Estimators->addColumn("Population");
    IndexType EtrialIndex = Estimators->addColumn("Etrial");
    Estimators->reportHeader(AppendRun);
    Estimators->reset();

    IndexType block = 0;
    RealType Eest = branchEngine->E_T;

    resetRun();

    for(int ip=0; ip<NumThreads; ip++) {
      char fname[16];
      sprintf(fname,"test.%i",ip);
      ofstream fout(fname);
    }

    for(int istep=0; istep<nSteps; istep++) {

      FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
#pragma omp parallel  
      {
        int ip = omp_get_thread_num();
        Movers[ip]->benchMark(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],ip);
      }
    }
    return true;
  }
  
  bool 
  DMCPbyPOpenMP::put(xmlNodePtr q){ 
    //nothing to do
    return true;
  }

}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
