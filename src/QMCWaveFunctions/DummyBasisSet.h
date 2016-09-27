//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
