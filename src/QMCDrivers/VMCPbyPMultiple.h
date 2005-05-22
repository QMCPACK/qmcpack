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
#ifndef OHMMS_QMC_VMC_PARTICLEBYPARTCLE_H
#define OHMMS_QMC_VMC_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCDriver.h" 
namespace ohmmsqmc {

  class MultipleEnergyEstimator;
  
  /** Implements the VMC algorithm using particle-by-particle move. */
  class VMCPbyPMultiple: public QMCDriver {
  public:
    /// Constructor.
    VMCPbyPMultiple(MCWalkerConfiguration& w, 
			  TrialWaveFunction& psi, 
			  QMCHamiltonian& h);
    bool run();
    bool put(xmlNodePtr cur);
 
  private:
    /// Copy Constructor (disabled)
    VMCPbyPMultiple(const VMCPbyPMultiple& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    VMCPbyPMultiple& operator=(const VMCPbyPMultiple&) { return *this;}

    int nPsi;
    typedef ParticleSet::ParticleGradient_t ParticleGradient_t;
    typedef ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;
    ParticleGradient_t dG;
    vector<ParticleGradient_t*> G;
    vector<ParticleLaplacian_t*> dL;
    vector<RealType> ratio, ratioij, UmbrellaWeight,sumratio;
    MultipleEnergyEstimator *multiEstimator;

    
    inline void resize(int n, int N){
      int m=n*(n-1); m=m/2;
      ratio.resize(n);
      UmbrellaWeight.resize(n);
      sumratio.resize(n);
      ratioij.resize(m);
      for(int i=0; i<n; i++){
	G.push_back(new ParticleGradient_t(N));
	dL.push_back(new ParticleLaplacian_t(N));
      }
    }
  };
  
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
