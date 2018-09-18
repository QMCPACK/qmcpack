//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file SPOSetProxy.cpp
 * @brief implements the member functions of SPOSetProxy
 */
#include "QMCWaveFunctions/Fermion/SPOSetProxyForMSD.h"
namespace qmcplusplus
{

template<Batching B>
void SPOSetProxyForMSD<B>::resetParameters(const opt_variables_type& optVariables)
{
  refPhi->resetParameters(optVariables);
}

template<Batching B>
void SPOSetProxyForMSD<B>::resetTargetParticleSet(ParticleSet& P)
{
  refPhi->resetTargetParticleSet(P);
}

template<Batching B>
void SPOSetProxyForMSD<B>::setOrbitalSetSize(int norbs)
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
  grad_grad_psiM.resize(OrbitalSetSize,norbs);
  grad_grad_grad_psiM.resize(OrbitalSetSize,norbs);
  grad_gradV.resize(norbs);
  grad_grad_gradV.resize(norbs);
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluateForPtclMove(const ParticleSet& P, int iat)
{
  refPhi->evaluate(P, iat, psiV);
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluateAllForPtclMove(const ParticleSet& P, int iat)
{
  refPhi->evaluate(P, iat, psiV, dpsiV, d2psiV);
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluateForWalkerMove(const ParticleSet& P, int first, int last)
{
  refPhi->evaluate_notranspose(P,first,last,psiM,dpsiM,d2psiM);
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluateForWalkerMoveWithHessian(const ParticleSet& P, int first, int last)
{
  refPhi->evaluate_notranspose(P,first,last,psiM,dpsiM,grad_grad_psiM);
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluateForWalkerMoveWithThirdDeriv(const ParticleSet& P, int first, int last)
{
  refPhi->evaluate_notranspose(P,first,last,psiM,dpsiM,grad_grad_psiM,grad_grad_grad_psiM);
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  int n=occup.cols();
  for(int i=0; i<n; i++)
  {
    int q = occup(workingSet,i);
    psi[i] = psiV[q];
  }
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluate(const ParticleSet& P, int iat,
                                 ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
{
  int n=occup.cols();
  for(int i=0; i<n; i++)
  {
    int q = occup(workingSet,i);
    psi[i] = psiV[q];
    dpsi[i] = dpsiV[q];
    d2psi[i] = d2psiV[q];
  }
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluate(const ParticleSet& P, int first, int last,
                                 ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  int n=occup.cols();
  for(int k=first,p=0; k<last; k++,p++)
  {
    for(int i=0; i<n; i++)
    {
      int q = occup(workingSet,i);
      logdet(i,p)=psiM(p,q);
      dlogdet(p,i)=dpsiM(p,q);
      d2logdet(p,i)=d2psiM(p,q);
    }
  }
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluate_notranspose(const ParticleSet& P, int first, int last,
    ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  int n=occup.cols();
  for(int k=first,p=0; k<last; k++,p++)
  {
    for(int i=0; i<n; i++)
    {
      int q = occup(workingSet,i);
      logdet(p,i)=psiM(p,q);
      dlogdet(p,i)=dpsiM(p,q);
      d2logdet(p,i)=d2psiM(p,q);
    }
  }
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluate(const ParticleSet& P, int first, int last,
                                 ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
{
  int n=occup.cols();
  for(int k=first,p=0; k<last; k++,p++)
  {
    for(int i=0; i<n; i++)
    {
      int q = occup(workingSet,i);
      logdet(i,p)=psiM(p,q);
      dlogdet(p,i)=dpsiM(p,q);
      grad_grad_logdet(p,i)=grad_grad_psiM(p,q);
    }
  }
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluate(const ParticleSet& P, int first, int last,
                                 ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet,
                                 GGGMatrix_t& grad_grad_grad_logdet)
{
  int n=occup.cols();
  for(int k=first,p=0; k<last; k++,p++)
  {
    for(int i=0; i<n; i++)
    {
      int q = occup(workingSet,i);
      logdet(i,p)=psiM(p,q);
      dlogdet(p,i)=dpsiM(p,q);
      grad_grad_logdet(p,i)=grad_grad_psiM(p,q);
      grad_grad_grad_logdet(p,i)=grad_grad_grad_psiM(p,q);
    }
  }
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluate_notranspose(const ParticleSet& P, int first, int last,
     ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
{
  int n=occup.cols();
  for(int k=first,p=0; k<last; k++,p++)
  {
    for(int i=0; i<n; i++)
    {
      int q = occup(workingSet,i);
      logdet(p,i)=psiM(p,q);
      dlogdet(p,i)=dpsiM(p,q);
      grad_grad_logdet(p,i)=grad_grad_psiM(p,q);
    }
  }
}

template<Batching B>
void SPOSetProxyForMSD<B>::evaluate_notranspose(const ParticleSet& P, int first, int last,
     ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet,
     GGGMatrix_t& grad_grad_grad_logdet)
{
  int n=occup.cols();
  for(int k=first,p=0; k<last; k++,p++)
  {
    for(int i=0; i<n; i++)
    {
      int q = occup(workingSet,i);
      logdet(p,i)=psiM(p,q);
      dlogdet(p,i)=dpsiM(p,q);
      grad_grad_logdet(p,i)=grad_grad_psiM(p,q);
      grad_grad_grad_logdet(p,i)=grad_grad_grad_psiM(p,q);
    }
  }
}

}
