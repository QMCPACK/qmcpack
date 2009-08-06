//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file SPOSetProxy.cpp
 * @brief implements the member functions of SPOSetProxy
 */
#include "QMCWaveFunctions/Fermion/SPOSetProxy.h"
namespace qmcplusplus {

  SPOSetProxy::SPOSetProxy(SPOSetBasePtr const& spos, int first, int last)
    : refPhi(spos)
  {
    Identity=true;
    classNamse="SPOSetProxy";
    OrbitalSetSize=last-first;
    BasisSetSize=last-first;
    setOrbitalSetSize(refPhi->getOrbitalSetSize());
  }

  void SPOSetProxy::resetParameters(const opt_variables_type& optVariables)
  {
    refPhi->resetParameters(optVariables);
  }

  void SPOSetProxy::resetTargetParticleSet(ParticleSet& P)
  {
    refPhi->resetTargetParticleSet(P);
  }

  void SPOSetProxy::setOrbitalSetSize(int norbs)
  {
    psiM.resize(norbs,OrbitalSetSize);
    dpsiM.resize(OrbitalSetSize,norbs);
    d2psiM.resize(OrbitalSetSize,norbs);
    psiV.resize(norbs);
    dpsiV.resize(norbs);
    d2psiV.resize(norbs);
  }

  void SPOSetProxy::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    Phi->evaluate(P, iat, psiV);
    std::copy(psiV.begin(),psiV.begin()+ObritalSetSize,psi.begin());
  }

  void SPOSetProxy::evaluate(const ParticleSet& P, int iat
      , ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    Phi->evaluate(P, iat, psiV,dpsiV,d2psiV);
    std::copy(psiV.begin(),psiV.begin()+ObritalSetSize,psi.begin());
    std::copy(dpsiV.begin(),dpsiV.begin()+ObritalSetSize,dpsi.begin());
    std::copy(d2psiV.begin(),d2psiV.begin()+ObritalSetSize,d2psi.begin());
  }

  void SPOSetProxy::evaluate(const ParticleSet& P, int first, int last
      , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    //evaluate all
    refPhi->evaluate(P,first,last,psiM,dpsiM,d2psiM);

    //copy the ground states
    std::copy(psiM.begin(),psiM.begin()+logdet.size(),logdet.begin());
    for(int i=0; i<OrbitalSetSize; ++i) std::copy(dpsiM[i],dpsiM[i]+OrbitalSetSize,dlogdet[i]);
    for(int i=0; i<OrbitalSetSize; ++i) std::copy(d2psiM[i],d2psiM[i]+OrbitalSetSize,d2logdet[i]);
  }

}
/***************************************************************************
 * $RCSfile$   $Author: kesler $
 * $Revision: 3535 $   $Date: 2009-02-10 13:04:12 -0600 (Tue, 10 Feb 2009) $
 * $Id: DeterminantTree.h 3535 2009-02-10 19:04:12Z kesler $ 
 ***************************************************************************/
