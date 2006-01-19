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
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"

namespace qmcplusplus { 

  /// Constructor.
  DMCPbyPOpenMP::DMCPbyPOpenMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, WaveFunctionFactory* psifac): 
    QMCDriver(w,psi,h),psiFactory(psifac){ 
    RootName = "dummy";
    QMCType ="dummy";
  }
  
  bool DMCPbyPOpenMP::run() { 
    
    m_oneover2tau = 0.5/Tau;
    m_sqrttau = sqrt(Tau);
    
    //cout << "Lattice of ParticleSet " << endl;
    //W.Lattice.print(cout);

    //property container to hold temporary properties, such as a local energy
    //MCWalkerConfiguration::PropertyContainer_t Properties;
    MCWalkerConfiguration::Walker_t& thisWalker(**W.begin());

    int np= omp_get_max_threads();
    cout << "The number of threads " << np << endl;
    vector<ParticleSet*> newPtclSets;
    vector<QMCHamiltonian*> Hdup;
    vector<TrialWaveFunction*> Psidup;

    //create arrays with 0 points
    newPtclSets.resize(np,0);
    Hdup.resize(np,0);
    Psidup.resize(np,0);
    psiFactory->setCloneSize(np);

    newPtclSets[0]=&W;
    Hdup[0]=&H;
    Psidup[0]=&Psi;


    char pnameIP[128];
#pragma omp parallel for
    for(int ip=0; ip<np; ip++) {
      if(newPtclSets[ip] == 0) {
        char pnameIP[128];
        sprintf(pnameIP,"%s.c%i",W.getName().c_str(),ip);
        newPtclSets[ip]=new ParticleSet(W);
        newPtclSets[ip]->setName(pnameIP);
        Hdup[ip]=new QMCHamiltonian(H);
        Hdup[ip]->resetTargetParticleSet(*newPtclSets[ip]);
        Psidup[ip]= psiFactory->cloneWaveFunction(newPtclSets[ip],ip);
      }
    }

    makeGaussRandom(deltaR); 

    drift = W.R;
    //new poosition
    //W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;

    vector<int> wPerNode(np+1);
    int n=W.getActiveWalkers()/np;
    wPerNode[0]=0;
    for(int ip=0; ip<np; ip++) {
      wPerNode[ip+1]=wPerNode[ip]+n;
    }
//#pragma omp parallel for
//    for(int ip=0; ip<np; ip++) {
    for(int istep=0; istep<nSteps; istep++) {
#pragma omp parallel  
    {
      int ip = omp_get_thread_num();
      char fname[16];
      sprintf(fname,"test%i.%i",ip,istep);
      ofstream fout(fname);
      //for(int i=0; i<10000; i++) {
      MCWalkerConfiguration::iterator it(W.begin()+wPerNode[ip]); 
      MCWalkerConfiguration::iterator it_end(W.begin()+wPerNode[ip+1]); 
      int i=0;
      while(it != it_end) {
        RealType x=static_cast<RealType>(i%10);
        Walker_t& thisWalker(**it);
        newPtclSets[ip]->R = x*m_sqrttau*deltaR + thisWalker.R;
        //update the distance table associated with W
        newPtclSets[ip]->update();
        //evaluate wave function
        //update the properties: note that we are getting \f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
        //ValueType logpsi(Psi.evaluateLog(*newPtclSets[ip]));
        ValueType logpsi(Psidup[ip]->evaluateLog(*newPtclSets[ip]));
        RealType e = Hdup[ip]->evaluate(*newPtclSets[ip]);
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
