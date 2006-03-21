//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Simone Chiesa
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
#ifndef QMCPLUSPLUS_VMC_PARTICLEBYPARTCLE_WARP_H
#define QMCPLUSPLUS_VMC_PARTICLEBYPARTCLE_WARP_H
#include "QMCDrivers/QMCDriver.h" 
#include "SpaceWarp.h"
namespace qmcplusplus {

  class MultipleEnergyEstimator;
  
  /** @ingroup QMCDrivers MultiplePsi ParticleByParticle
   * @brief Implements the VMC algorithm
   */
  class VMCPbyPMultiWarp: public QMCDriver {
  public:
    /// Constructor.
    VMCPbyPMultiWarp(MCWalkerConfiguration& w, 
			  TrialWaveFunction& psi, 
			  QMCHamiltonian& h);
    ~VMCPbyPMultiWarp();

    bool run();
    bool put(xmlNodePtr cur);
 
  private:
    /// Copy Constructor (disabled)
    VMCPbyPMultiWarp(const VMCPbyPMultiWarp& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    VMCPbyPMultiWarp& operator=(const VMCPbyPMultiWarp&) { return *this;}

    int nPsi,nptcl,JACOBIAN;
    typedef ParticleSet::ParticleGradient_t ParticleGradient_t;
    typedef ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;
    ParticleGradient_t dG;
    vector<ParticleGradient_t*> G;
    vector<ParticleLaplacian_t*> dL;
    vector<RealType> ratio, ratioij, logpsi2, UmbrellaWeight,sumratio,invsumratio;
    vector<ParticleSet*> WW;
    MultipleEnergyEstimator *multiEstimator;

    SpaceWarp PtclWarp;
    
    inline void resize(int ncopy, int nptcls){
      int m=ncopy*(ncopy-1)/2; 
      ratio.resize(ncopy);
      logpsi2.resize(ncopy);
      UmbrellaWeight.resize(ncopy);
      invsumratio.resize(ncopy);
      sumratio.resize(ncopy);
      ratioij.resize(m);
      for(int i=0; i<ncopy; i++){
	G.push_back(new ParticleGradient_t(nptcls));
	dL.push_back(new ParticleLaplacian_t(nptcls));
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
