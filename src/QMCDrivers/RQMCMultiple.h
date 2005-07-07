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
#ifndef OHMMS_QMC_REPMULTIPLE_H
#define OHMMS_QMC_REPMULTIPLE_H

#include "QMCDrivers/QMCDriver.h" 
#include "OhmmsPETE/OhmmsVector.h" 

namespace ohmmsqmc {

  class PolymerChain;

  /** Implements the RMC algorithm. */
  class RQMCMultiple: public QMCDriver {

  public:

    /// Constructor.
    RQMCMultiple(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

    /// Destructor
    ~RQMCMultiple();

    bool run();
    bool put(xmlNodePtr q);

  protected:

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

    ///The length of polymers
    int ReptileLength;

    ///the number of the beads that will be cut
    int  NumCuts;

    ///the number of turns per block
    int NumTurns;

    ///the number of H/Psi pairs
    int nPsi;

    int ForwardKineticAction;
    int BackwardKineticAction;

    ///array of PolymerChains
    PolymerChain* Reptile;

    ///move polymers
    void moveReptile();

    ///initialize polymers
    void initReptile();

    ///Working arrays
    Vector<RealType> SumRatioAction,LogRatioActionIJ,sumratio,WReptile,logpsi;
    Vector<int> TotalSign,SignWeight,beadSignWeight;

    ParticleSet::ParticlePos_t DiffusionDrift;
    ParticleSet::ParticleGradient_t LocalDrift;

    void resizeArrays(int n);

  private:

    /// Copy Constructor (disabled)
    RQMCMultiple(const RQMCMultiple& a): QMCDriver(a) { }

    /// Copy operator (disabled).
    RQMCMultiple& operator=(const RQMCMultiple&) { return *this;}


  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
