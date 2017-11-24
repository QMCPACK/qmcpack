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
    
    
/** @file SPOSetProxy.h
 * @brief declare a proxy class to a SPOSetBase for multi determinants
 */
#ifndef QMCPLUSPLUS_SPOSETPROXY_FORMSD_H
#define QMCPLUSPLUS_SPOSETPROXY_FORMSD_H
#include "QMCWaveFunctions/SPOSetBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"
namespace qmcplusplus
{

/** proxy SPOSetBase
 *
 * This class owns a SPOSetBase for all the states to be evaluated
 * and will be owned by a DiracDeterminantBase object.
 */
struct SPOSetProxyForMSD: public SPOSetBase
{

  ///pointer to the SPOSet which evaluate the single-particle states
  SPOSetBasePtr refPhi;
  ///container for the values
  ValueMatrix_t psiM;
  ///container for the gradients
  GradMatrix_t dpsiM;
  ///container for the laplacians
  ValueMatrix_t d2psiM;
  ///container for the hessians
  HessMatrix_t grad_grad_psiM;
  ///container for the hessian
  HessVector_t grad_gradV;
  ///container for the GGGtypes
  GGGMatrix_t  grad_grad_grad_psiM;
  ///container for the GGGtypes
  GGGVector_t  grad_grad_gradV;
  ///contatiner for the values for a particle
  ValueVector_t psiV;
  ///contatiner for the gradients for a particle
  GradVector_t dpsiV;
  ///contatiner for the laplacians for a particle
  ValueVector_t d2psiV;

  Matrix<int> occup;

  int workingSet;

  /** constructor
   * @param spos a SPOSet
   * @param first the first particle index
   * @param last the last particle index
   */
  SPOSetProxyForMSD(SPOSetBasePtr const& spos, int first, int last);
  void resetParameters(const opt_variables_type& optVariables);
  void resetTargetParticleSet(ParticleSet& P);
  void setOrbitalSetSize(int norbs);
  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat
                , ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
  inline void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& d2psi)
  {
    APP_ABORT("Need specialization of evaluate(HessVector_t&) for grad_grad_grad_logdet. \n");
  }
  void evaluate(const ParticleSet& P, int first, int last
                , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
  void evaluate(const ParticleSet& P, int first, int last
                , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);

  void evaluate(const ParticleSet& P, int first, int last
                , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet);

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet);

  void prepareFor(int n)
  {
    workingSet = n;
  }

  void evaluateForWalkerMove(const ParticleSet& P, int first, int last);
  void evaluateForWalkerMoveWithHessian(const ParticleSet& P, int first, int last);
  void evaluateForWalkerMoveWithThirdDeriv(const ParticleSet& P, int first, int last);
  void evaluateForPtclMove(const ParticleSet& P, int iat);
  void evaluateAllForPtclMove(const ParticleSet& P, int iat);

};
}
#endif
