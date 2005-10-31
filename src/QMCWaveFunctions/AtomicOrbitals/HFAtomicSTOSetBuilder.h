//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_ATOMIC_HARTREEFOCK_STO_IOXML_H
#define QMCPLUSPLUS_ATOMIC_HARTREEFOCK_STO_IOXML_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSet.h"

namespace qmcplusplus {

  /** OrbitalBuilder with HFAtomicSTOSet
   */
  class HFAtomicSTOSetBuilder: public OrbitalBuilderBase {

  private:

    typedef HFAtomicSTOSet::RadialOrbital_t RadialOrbital_t;
    typedef HFAtomicSTOSet::SPO_t           SPO_t;
    DistanceTableData* d_table;

    /// Maximum Angular Momentum
    int Lmax;

    /// mapping function for Rnl[ RnlID[name] ]
    map<string,int> RnlID;         

    ///temporary storage of radial functions
    vector<RadialOrbital_t*> Rnl;            

    ///single-particle wave functions
    map<string, SPO_t*>   OrbSet;

    bool getBasis(xmlNodePtr cur);
    HFAtomicSTOSet* getOrbital(xmlNodePtr cur);

  public:

    HFAtomicSTOSetBuilder(ParticleSet& els, TrialWaveFunction& psi, ParticleSet& ions);
    bool put(xmlNodePtr);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
