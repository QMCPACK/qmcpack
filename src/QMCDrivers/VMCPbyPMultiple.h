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
#ifndef OHMMS_QMC_VMC_PARTICLEBYPARTCLE_MULTIPLE_H
#define OHMMS_QMC_VMC_PARTICLEBYPARTCLE_MULTIPLE_H
#include "QMCDrivers/QMCDriver.h" 
namespace ohmmsqmc {

  class MultipleEnergyEstimator;
  
  /** @ingroup QMCDrivers MultiplePsi ParticleByParticle
   * @brief Implements the VMC algorithm
   */
  class VMCPbyPMultiple: public QMCDriver {
  public:
    /// Constructor.
    VMCPbyPMultiple(MCWalkerConfiguration& w, 
			  TrialWaveFunction& psi, 
			  QMCHamiltonian& h);
    ~VMCPbyPMultiple();

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
    vector<RealType> ratio, ratioij, logpsi2, UmbrellaWeight,sumratio,invsumratio;
    MultipleEnergyEstimator *multiEstimator;

    
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
