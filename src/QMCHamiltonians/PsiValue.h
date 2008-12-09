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
#ifndef QMCPLUSPLUS_PSIVALUE_H
#define QMCPLUSPLUS_PSIVALUE_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus {

  struct PsiValue: public QMCHamiltonianBase {

    PsiValue( ){ 
      UpdateMode.set(OPTIMIZABLE,1);
    }
    ///destructor
    ~PsiValue() { }

    void resetTargetParticleSet(ParticleSet& P) { }

    inline Return_t 
    evaluate(ParticleSet& P) {
      Value = std::exp(-2.0*tWalker->Properties(LOGPSI));
      return Value;
    }

    inline Return_t 
    evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    inline Return_t 
    registerData(ParticleSet& P, BufferType& buffer) 
    {
      Value = std::exp(tWalker->Properties(LOGPSI));
      buffer.add(Value);
      return Value;
    }

    inline Return_t 
    updateBuffer(ParticleSet& P, BufferType& buffer) 
    {
      Value = std::exp(tWalker->Properties(LOGPSI));
      buffer.put(Value);
      return Value;
    }

    inline void copyFromBuffer(ParticleSet& P, BufferType& buffer)
    {
      buffer.get(Value);
    }

    inline void copyToBuffer(ParticleSet& P, BufferType& buffer)
    {
      buffer.put(Value);
    }

    inline Return_t 
    evaluatePbyP(ParticleSet& P, int active)
    {
      Value = std::exp(tWalker->Properties(LOGPSI));
      return NewValue;
    }
 
    bool put(xmlNodePtr) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "PsiValue";
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return new PsiValue(*this);
    }

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3111 $   $Date: 2008-09-14 13:20:23 -0500 (Sun, 14 Sep 2008) $
 * $Id: PsiValue.h 3111 2008-09-14 18:20:23Z jnkim $ 
 ***************************************************************************/

