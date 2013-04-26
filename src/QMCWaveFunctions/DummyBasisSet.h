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

namespace qmcplusplus
{

//  struct BasisSetBase {
//    virtual void resize(int ntpcl) = 0;
//    virtual void resetTargetParticleSet(ParticleSet& P)=0;
//    virtual void evaluate(const ParticleSet& P)=0;
//    virtual void evaluate(const ParticleSet& P, int iat)=0;
//    virtual void evaluateAll(const ParticleSet& P, int iat)=0;
//  };
//
//  template<class BS>
//  struct BasisSetProxy: public BasisSetBase {
//    BS* basisRef;
//    BasisSetProxy(BS* basis): basisRef(basis){}
//
//    void resize(int nptcl) { basisRef->resize(nptcl);}
//
//    void evaluate(const ParticleSet& P) {
//      basisRef->evaluate(P);
//    }
//    void evaluate(const ParticleSet& P, int iat) {
//      basisRef->evaluate(P,iat);
//    }
//    void evaluateAll(const ParticleSet& P, int iat) {
//      basisRef->evaluateAll(P,iat);
//    }
//  };
//
/** DummyBasisSet class which implements functions
 *
 * Any S(ingle)P(article)O(rbital)Set can contain
 * typedef DummyBasisSet BasisSet_t
 * so that SlaterDeterminant<SPOSet> can call the dummy function.
 *
 */
struct DummyBasisSet
{
  inline void resetTargetParticleSet(ParticleSet& P) { }
  inline void evaluate(ParticleSet& P) { }
};


}
#endif
