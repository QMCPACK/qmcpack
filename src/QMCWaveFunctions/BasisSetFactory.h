//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_BASISSETFACTORY_H
#define QMCPLUSPLUS_BASISSETFACTORY_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/BasisSetBase.h"

namespace qmcplusplus {

  /** derived class from OrbitalBuilderBase
   */
  class BasisSetFactory: public OrbitalBuilderBase {

  public:

    /** constructor
     * \param els reference to the electrons
     * \param psi reference to the wavefunction
     * \param ions reference to the ions
     */
    BasisSetFactory(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets);

    bool put(xmlNodePtr cur);

    void createBasisSet(xmlNodePtr cur, xmlNodePtr rootNode=NULL);

    SPOSetBase* createSPOSet(xmlNodePtr cur);

  private:
    ///set of basis set: potential static data
    vector<BasisSetBuilder*> basisBuilder;
    ///reference to the particle pool
    PtclPoolType& ptclPool;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
