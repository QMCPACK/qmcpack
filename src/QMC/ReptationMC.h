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
#ifndef OHMMS_QMC_REPTATION_H
#define OHMMS_QMC_REPTATION_H

#include "QMC/QMCDriver.h" 
#include <deque>
namespace ohmmsqmc {


  class PolymerChain;

  /** Implements the RMC algorithm. */
  class ReptationMC: public QMCDriver {

  public:

    /// Constructor.
    ReptationMC(MCWalkerConfiguration& w, 
	       TrialWaveFunction& psi, 
	       QMCHamiltonian& h, 
	       xmlNodePtr q);

    /// Destructor
    ~ReptationMC();

    bool run();
    bool put(xmlNodePtr q);

  protected:

    typedef MCWalkerConfiguration::Walker_t Walker_t;

    ///boolean for using bounce algorithm. true, if bounce algorithm of D. Ceperley
    bool UseBounce;

    /** boolean for initialization
     *
     *\if true,
     *use clones for a chain.
     *\else
     *use drift-diffusion to form a chain
     *\endif
     */
    bool ClonePolymer;

    ///The length of polymers
    int PolymerLength;

    ///the number of the beads that will be cut
    int  NumCuts;

    ///the number of turns per block
    int NumTurns;

    ///temporary storage for random displacement
    ParticleSet::ParticlePos_t deltaR;

    ///array of PolymerChains
    std::vector<PolymerChain*> Polymers;

    ///move polymers
    void movePolymers();

    ///initialize polymers
    void initPolymers();
  private:

    /// Copy Constructor (disabled)
    ReptationMC(const ReptationMC& a): QMCDriver(a) { }

    /// Copy operator (disabled).
    ReptationMC& operator=(const ReptationMC&) { return *this;}


  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
