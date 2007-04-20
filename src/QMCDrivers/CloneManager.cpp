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
#include "QMCDrivers/CloneManager.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus { 

  //initialization of the static wClones
  vector<ParticleSet*> CloneManager::wClones;
  //initialization of the static psiClones
  vector<TrialWaveFunction*> CloneManager::psiClones;
  //initialization of the static hClones
  vector<QMCHamiltonian*> CloneManager::hClones;

  /// Constructor.
  CloneManager::CloneManager(HamiltonianPool& hpool): cloneEngine(hpool)
  {
    NumThreads=omp_get_max_threads();
    wPerNode.resize(NumThreads+1,0);
  }

  ///clenup non-static data members
  CloneManager::~CloneManager()
  {
    delete_iter(Rng.begin(),Rng.end());
    delete_iter(Movers.begin(),Movers.end());
    delete_iter(branchClones.begin(),branchClones.end());
    delete_iter(estimatorClones.begin(),estimatorClones.end());
  }

  void CloneManager::makeClones(MCWalkerConfiguration& w, 
      TrialWaveFunction& psi, QMCHamiltonian& ham)
  {

    if(wClones.size()) {
      app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << endl;
      return;
    }

    app_log() << "Number of threads = " << NumThreads << endl;
    wClones.resize(NumThreads,0);
    psiClones.resize(NumThreads,0);
    hClones.resize(NumThreads,0);

    wClones[0]=&w;
    psiClones[0]=&psi;
    hClones[0]=&ham;

    cloneEngine.clone(w,psi,ham,wClones,psiClones,hClones);
  }
}

/***************************************************************************
 * $RCSfile: DMCPbyPOMP.cpp,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/10/18 17:23:35 $
 * $Id: DMCPbyPOMP.cpp,v 1.5 2006/10/18 17:23:35 jnkim Exp $ 
 ***************************************************************************/
