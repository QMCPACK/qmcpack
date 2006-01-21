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
    wClones[0]=&W;
    psiClones[0]=&Psi;
    hClones[0]=&H;
    if(NumThreads == 1) {
      WARNMSG("Using a single thread with DMCPbyPOpenMP.")
      return;
    }
    hpool.clone(W,Psi,H,wClones,psiClones,hClones);
  }
  
  bool DMCPbyPOpenMP::run() { 
    
    m_oneover2tau = 0.5/Tau;
    m_sqrttau = sqrt(Tau);
    
    //cout << "Lattice of ParticleSet " << endl;
    //W.Lattice.print(cout);

    //property container to hold temporary properties, such as a local energy
    //MCWalkerConfiguration::PropertyContainer_t Properties;
    MCWalkerConfiguration::Walker_t& thisWalker(**W.begin());

    makeGaussRandom(deltaR); 

    drift = W.R;
    //new poosition
    //W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;

    int n=W.getActiveWalkers()/NumThreads;
    for(int ip=0; ip<NumThreads; ip++) {
      wPerNode[ip+1]=wPerNode[ip]+n;
    }

    
    for(int ip=0; ip<NumThreads; ip++) {
        char fname[16];
        sprintf(fname,"test.%i",ip);
        ofstream fout(fname);
    }
//#pragma omp parallel for
//    for(int ip=0; ip<np; ip++) {

    //TESTING
    for(int istep=0; istep<nSteps; istep++) {
#pragma omp parallel  
      {
        int ip = omp_get_thread_num();
        char fname[16];
        sprintf(fname,"test.%i",ip);
        ofstream fout(fname,ios::app);
        //for(int i=0; i<10000; i++) {
        MCWalkerConfiguration::iterator it(W.begin()+wPerNode[ip]); 
        MCWalkerConfiguration::iterator it_end(W.begin()+wPerNode[ip+1]); 
        int i=0;
        while(it != it_end) {
          RealType x=static_cast<RealType>(i%10);
          Walker_t& thisWalker(**it);
          wClones[ip]->R = x*m_sqrttau*deltaR + thisWalker.R;
          //update the distance table associated with W
          wClones[ip]->update();
          //evaluate wave function
          //update the properties: note that we are getting \f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
          //ValueType logpsi(Psi.evaluateLog(*newPtclSets[ip]));
          ValueType logpsi(psiClones[ip]->evaluateLog(*wClones[ip]));
          RealType e = hClones[ip]->evaluate(*wClones[ip]);
          fout << logpsi << " " << e << endl;
          ++i; ++it;
        }
        //}
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
