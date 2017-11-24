//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SPOSetProxy.cpp
 * @brief implements the member functions of SPOSetProxy
 */
#include "QMCWaveFunctions/Fermion/SPOSetProxy.h"
namespace qmcplusplus
{

SPOSetProxy::SPOSetProxy(SPOSetBasePtr const& spos, int first, int last)
  : refPhi(spos)
{
  Identity=true;
  className="SPOSetProxy";
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
  //psiM.resize(norbs,OrbitalSetSize);
  //dpsiM.resize(norbs,OrbitalSetSize);
  //d2psiM.resize(norbs,OrbitalSetSize);
  psiM.resize(OrbitalSetSize,norbs);
  dpsiM.resize(OrbitalSetSize,norbs);
  d2psiM.resize(OrbitalSetSize,norbs);
  psiV.resize(norbs);
  dpsiV.resize(norbs);
  d2psiV.resize(norbs);
}

void SPOSetProxy::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  refPhi->evaluate(P, iat, psiV);
  std::copy(psiV.begin(),psiV.begin()+OrbitalSetSize,psi.begin());
// mmorales: needed for MultiSlaterDeterminant moves: put an if statement??
  std::copy(psiV.begin(),psiV.end(),psiM[iat]);
}

void SPOSetProxy::evaluate(const ParticleSet& P, int iat,
                           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
{
  refPhi->evaluate(P, iat, psiV, dpsiV, d2psiV);
  std::copy(psiV.begin(),psiV.begin()+OrbitalSetSize,psi.begin());
  std::copy(dpsiV.begin(),dpsiV.begin()+OrbitalSetSize,dpsi.begin());
  std::copy(d2psiV.begin(),d2psiV.begin()+OrbitalSetSize,d2psi.begin());
// mmorales: needed for MultiSlaterDeterminant moves: put an if statement??
  std::copy(psiV.begin(),psiV.end(),psiM[iat]);
  std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[iat]);
  std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[iat]);
}

void SPOSetProxy::evaluate(const ParticleSet& P, int first, int last
                           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  //evaluate all using notranspose
  refPhi->evaluate_notranspose(P,first,last,psiM,dpsiM,d2psiM);
  //transpose only the ground state
  for(int i=0; i<OrbitalSetSize; ++i)
    for(int j=0; j<OrbitalSetSize; ++j)
      logdet(i,j)=psiM(j,i);
  //copy the ground states
  for(int i=0; i<OrbitalSetSize; ++i)
    std::copy(dpsiM[i],dpsiM[i]+OrbitalSetSize,dlogdet[i]);
  for(int i=0; i<OrbitalSetSize; ++i)
    std::copy(d2psiM[i],d2psiM[i]+OrbitalSetSize,d2logdet[i]);
}

void SPOSetProxy::evaluate_notranspose(const ParticleSet& P, int first, int last
                                       , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  //evaluate all
  refPhi->evaluate_notranspose(P,first,last,psiM,dpsiM,d2psiM);
  for(int i=0; i<OrbitalSetSize; ++i)
    for(int j=0; j<OrbitalSetSize; ++j)
      logdet(i,j)=psiM(i,j);
  for(int i=0; i<OrbitalSetSize; ++i)
    std::copy(dpsiM[i],dpsiM[i]+OrbitalSetSize,dlogdet[i]);
  for(int i=0; i<OrbitalSetSize; ++i)
    std::copy(d2psiM[i],d2psiM[i]+OrbitalSetSize,d2logdet[i]);
}

void SPOSetProxy::evaluate(const ParticleSet& P, int first, int last
                           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
{
  APP_ABORT("SPOSetProxy::evaluate_notranspose need specialization for GGGMatrix.\n");
}

void SPOSetProxy::evaluate(const ParticleSet& P, int first, int last
                           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet
                           , GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("SPOSetProxy::evaluate_notranspose need specialization for GGGMatrix.\n");
}

void SPOSetProxy::evaluate_notranspose(const ParticleSet& P, int first, int last
                                       , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
{
  APP_ABORT("SPOSetProxy::evaluate_notranspose need specialization for GGGMatrix.\n");
}

void SPOSetProxy::evaluate_notranspose(const ParticleSet& P, int first, int last
                                       , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet
                                       , GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("SPOSetProxy::evaluate_notranspose need specialization for GGGMatrix.\n");
}

}
