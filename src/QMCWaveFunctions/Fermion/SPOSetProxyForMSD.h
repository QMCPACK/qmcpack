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
 * @brief declare a proxy class to a SPOSet for multi determinants
 */
#ifndef QMCPLUSPLUS_SPOSETPROXY_FORMSD_H
#define QMCPLUSPLUS_SPOSETPROXY_FORMSD_H
#include "QMCWaveFunctions/SPOSet.h"
#include "OhmmsPETE/OhmmsMatrix.h"
namespace qmcplusplus
{
/** proxy SPOSet
 *
 * This class owns a SPOSet for all the states to be evaluated
 * and will be owned by a DiracDeterminant object.
 */
struct SPOSetProxyForMSD : public SPOSet
{
  ///pointer to the SPOSet which evaluate the single-particle states
  std::unique_ptr<SPOSet> refPhi;
  ///container for the values
  ValueMatrix psiM;
  ///container for the gradients
  GradMatrix dpsiM;
  ///container for the laplacians
  ValueMatrix d2psiM;
  ///container for the hessians
  HessMatrix grad_grad_psiM;
  ///container for the hessian
  HessVector grad_gradV;
  ///container for the GGGtypes
  GGGMatrix grad_grad_grad_psiM;
  ///container for the GGGtypes
  GGGVector grad_grad_gradV;
  ///contatiner for the values for a particle
  ValueVector psiV;
  ///contatiner for the gradients for a particle
  GradVector dpsiV;
  ///contatiner for the laplacians for a particle
  ValueVector d2psiV;

  Matrix<int> occup;

  int workingSet;

  /** constructor
   * @param spos a SPOSet
   * @param first the first particle index
   * @param last the last particle index
   */
  SPOSetProxyForMSD(std::unique_ptr<SPOSet>&& spos, int first, int last);
  void resetParameters(const opt_variables_type& optVariables) override;
  void setOrbitalSetSize(int norbs) override;
  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;
  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;
  void evaluate(const ParticleSet& P,
                int first,
                int last,
                ValueMatrix& logdet,
                GradMatrix& dlogdet,
                ValueMatrix& d2logdet);
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;
  void evaluate(const ParticleSet& P,
                int first,
                int last,
                ValueMatrix& logdet,
                GradMatrix& dlogdet,
                HessMatrix& grad_grad_logdet);

  void evaluate(const ParticleSet& P,
                int first,
                int last,
                ValueMatrix& logdet,
                GradMatrix& dlogdet,
                HessMatrix& grad_grad_logdet,
                GGGMatrix& grad_grad_grad_logdet);

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet,
                            GGGMatrix& grad_grad_grad_logdet) override;

  void prepareFor(int n) { workingSet = n; }

  void evaluateForWalkerMove(const ParticleSet& P, int first, int last);
  void evaluateForWalkerMoveWithHessian(const ParticleSet& P, int first, int last);
  void evaluateForWalkerMoveWithThirdDeriv(const ParticleSet& P, int first, int last);
  void evaluateForPtclMove(const ParticleSet& P, int iat);
  void evaluateAllForPtclMove(const ParticleSet& P, int iat);
};
} // namespace qmcplusplus
#endif
