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
#ifndef OHMMS_QMC_VMCMultiple_H
#define OHMMS_QMC_VMCMultiple_H
#include "QMC/QMCDriver.h" 
namespace ohmmsqmc {

  /** Implements the VMCMultiple algorithm. 
   */
  class VMCMultiple: public QMCDriver {
  public:
    /// Constructor.
    VMCMultiple(MCWalkerConfiguration& w, 
	TrialWaveFunction& psi, // SIMONE - Is this needed if we 
	QMCHamiltonian& h,      //          are goint to use Psi1 ?
	xmlNodePtr q);

    void advanceWalkerByWalker();
    bool run();
    bool put(xmlNodePtr cur);
 
  private:
    /// Copy Constructor (disabled)
    VMCMultiple(const VMCMultiple& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    VMCMultiple& operator=(const VMCMultiple&) { return *this;}

    ///temporary storage for drift
    //ParticleSet::ParticlePos_t drift;		//SIMONE
    typedef ParticleSet::ParticlePos_t       ParticlePos_t;	//SIMONE
    typedef ParticleSet::ParticleGradient_t  ParticleGradient_t;	//SIMONE
    typedef ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;	//SIMONE
    ParticlePos_t drift;				//SIMONE

    ///temporary storage for random displacement
    ParticlePos_t deltaR;

    ///temporary storage
    int nPsi;				//SIMONE
    vector<IndexType>IndexPsi;		//SIMONE
    vector<RealType> logpsi;		//SIMONE
    vector<ParticleGradient_t>  dgrad;		//SIMONE
    vector<RealType> sumratio;		//SIMONE
    vector<ParticleLaplacian_t> lap;		//SIMONE
    vector<RealType> invsumratio;	//SIMONE
    vector<RealType> totinvsumratio;	//SIMONE
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
