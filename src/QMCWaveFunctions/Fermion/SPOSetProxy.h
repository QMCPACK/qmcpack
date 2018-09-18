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
    
    
/** @file SPOSetProxy.h
 * @brief declare a proxy class to a SPOSet for multi determinants
 */
#ifndef QMCPLUSPLUS_SPOSETPROXY_H
#define QMCPLUSPLUS_SPOSETPROXY_H
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/SPOSetBatched.h"

namespace qmcplusplus
{

/** proxy SPOSet
 *
 * This class owns a SPOSet for all the states to be evaluated
 * and will be owned by a DiracDeterminantBase object.
 */
template<Batching B>
struct SPOSetProxy : public SPOSet<B>
{
  using SPOSetPtr = SPOSet<B>*;
  ///pointer to the SPOSet which evaluate the single-particle states
  SPOSetPtr refPhi;
  using SSTA = SPOSetTypeAliases;
  using ValueVector_t = SSTA::ValueVector_t;  
  using ValueMatrix_t = SSTA::ValueMatrix_t;
  using GradVector_t = SSTA::GradVector_t;
  using GradMatrix_t = SSTA::GradMatrix_t;
  using HessVector_t = SSTA::HessVector_t;
  using HessMatrix_t = SSTA::HessMatrix_t;
  using HessType = SSTA::HessType;
  using HessArray_t = SSTA::HessArray_t;
  using GGGMatrix_t = SSTA::GGGMatrix_t;
    
  ///container for the values
  ValueMatrix_t psiM;
  ///container for the gradients
  GradMatrix_t dpsiM;
  ///container for the laplacians
  ValueMatrix_t d2psiM;
  ///contatiner for the values for a particle
  ValueVector_t psiV;
  ///contatiner for the gradients for a particle
  GradVector_t dpsiV;
  ///contatiner for the laplacians for a particle
  ValueVector_t d2psiV;

  /** constructor
   * @param spos a SPOSet<Batching::SINGLE>
   * @param first the first particle index
   * @param last the last particle index
   */
SPOSetProxy<B>::SPOSetProxy(SPOSetPtr const& spos, int first, int last)
  : refPhi(spos)
{
  className="SPOSetProxy";
  OrbitalSetSize=last-first;
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


  /* SPOSetProxy(SPOSetPtr const& spos, int first, int last); */
  /* void resetParameters(const opt_variables_type& optVariables); */
  /* void resetTargetParticleSet(ParticleSet& P); */
  /* void setOrbitalSetSize(int norbs); */
  /* void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi); */
  /* void evaluate(const ParticleSet& P, int iat, */
  /*               ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi); */
  inline void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& d2psi)
  {
    APP_ABORT("Need specialization of evaluate(HessVector_t&) for grad_grad_grad_logdet. \n");
  }
  /* void evaluate(const ParticleSet& P, int first, int last */
  /*               , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet); */
  /* void evaluate_notranspose(const ParticleSet& P, int first, int last */
  /*                           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet); */
  /* void evaluate(const ParticleSet& P, int first, int last */
  /*               , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet); */

  /* void evaluate(const ParticleSet& P, int first, int last */
  /*               , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet); */

  /* void evaluate_notranspose(const ParticleSet& P, int first, int last */
  /*                           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet); */

  /* void evaluate_notranspose(const ParticleSet& P, int first, int last */
  /*                           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet); */

};
}
#endif
