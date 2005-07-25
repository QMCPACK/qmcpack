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

  class Bead;
  
  class MultiChain;

  /** @ingroup QMCDrivers MultiplePsi
   * @brief Implements the RMC algorithm for energy differences
   */
  class RQMCMultiple: public QMCDriver {

  public:

    /// Constructor.
    RQMCMultiple(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

    /// Destructor
    ~RQMCMultiple();

    bool run();
    bool put(xmlNodePtr q);

  protected:

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

    ///
    int MinusDirection,PlusDirection,Directionless;
    ///
    int forward,backward,ianchor,itail,inext;

    ///the number of turns per block
    int NumTurns;

    int nptcl;

    ///the number of H/Psi pairs
    int nPsi;

    ///The Reptile: a chain of beads
    MultiChain* Reptile;

    ///The new bead
   Bead *NewBead;

    ///move polymers
    void moveReptile();

    ///initialize polymers
    void initReptile();

    ///Working arrays
    Vector<RealType> SumRatioAction,LogRatioActionIJ,sumratio,WReptile,logpsi;

    Vector<int>NewTotalSign,WeightSign,RefSign;
    
    void resizeArrays(int n);

  private:

    /// Copy Constructor (disabled)
    RQMCMultiple(const RQMCMultiple& a): QMCDriver(a) { }

    /// Copy operator (disabled).
    RQMCMultiple& operator=(const RQMCMultiple&) { return *this;}

    double NewLogRatioAction,NewSumRatio,accept;

    ParticleSet::ParticlePos_t gRand;

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
