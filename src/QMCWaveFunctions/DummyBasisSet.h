//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_DUMMYBASISSET_H
#define QMCPLUSPLUS_DUMMYBASISSET_H

#include "Particle/ParticleSet.h"

namespace qmcplusplus {
  /** DummyBasisSet class which implements functions
   *
   * Any S(ingle)P(article)O(rbital)Set can contain
   * typedef DummyBasisSet BasisSet_t
   * so that SlaterDeterminant<SPOSet> can call the dummy function.
   *
   */
  struct DummyBasisSet {
    inline void resetTargetParticleSet(ParticleSet& P) { }
    inline void evaluate(ParticleSet& P) { }
  };

  struct BasisSetBase {
    virtual void resize(int n) = 0;
  };

  template<class BS>
  struct BasisSetProxy: public BasisSetBase {
    BS* basisRef;
    BasisSetProxy(BS* basis): basisRef(basis){}
    void resize(int n) { basisRef->resize(n);}
  };
}
#endif
