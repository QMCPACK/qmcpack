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
namespace qmcplusplus
{
/** proxy SPOSet
 *
 * This class owns a SPOSet for all the states to be evaluated
 * and will be owned by a DiracDeterminant object.
 */
struct SPOSetProxy : public SPOSet
{
  ///pointer to the SPOSet which evaluate the single-particle states
  SPOSetPtr refPhi;
  ///container for the values
  ValueMatrix psiM;
  ///container for the gradients
  GradMatrix dpsiM;
  ///container for the laplacians
  ValueMatrix d2psiM;
  ///contatiner for the values for a particle
  ValueVector psiV;
  ///contatiner for the gradients for a particle
  GradVector dpsiV;
  ///contatiner for the laplacians for a particle
  ValueVector d2psiV;

  /** constructor
   * @param spos a SPOSet
   * @param first the first particle index
   * @param last the last particle index
   */
  SPOSetProxy(SPOSetPtr const& spos, int first, int last);
  void resetParameters(const opt_variables_type& optVariables) override;
  void setOrbitalSetSize(int norbs) override;
  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;
  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;
  inline void evaluateVGH(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, HessVector& d2psi) override
  {
    APP_ABORT("Need specialization of evaluate(HessVector&) for grad_grad_grad_logdet. \n");
  }
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
};
} // namespace qmcplusplus
#endif
