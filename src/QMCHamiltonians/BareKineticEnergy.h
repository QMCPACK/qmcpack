//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_BAREKINETICENERGY_H
#define QMCPLUSPLUS_BAREKINETICENERGY_H

#include "QMCHamiltonians/OperatorBase.h"
#include <ResourceCollection.h>
#include <ResourceHandle.h>

namespace qmcplusplus
{

/** @ingroup hamiltonian
  @brief Evaluate the kinetic energy with a single mass

 *The unit of the mass is AU, i.e., the electron mass \f$ m_e = 1 \f$.
 * To evaluate the Bare Kinetic part of the local energy
 \f$E_L({\bf R}) = \Psi^{-1}({\bf R})\hat{H}\Psi({\bf R}),\f$
 it is useful to use the following trick
 \f{eqnarray*}
 \nabla^2\Psi({\bf R}) &=& \nabla^2(\exp(\ln \Psi({\bf R})))\\
 &=&\nabla\cdot(\nabla\exp(\ln \Psi({\bf R}))) \\
 &=&\nabla\cdot(\nabla\ln \Psi({\bf R}))\exp(\ln \Psi({\bf R}))\\
 -\frac{1}{2}\frac{\nabla^2\Psi({\bf R})}{\Psi({\bf R})} &=&
 -\frac{1}{2}\nabla^2\ln \Psi({\bf R})
 -\frac{1}{2}(\nabla\ln \Psi({\bf R}))^2
 \f}
 */
class BareKineticEnergy : public OperatorBase
{
public:
  /** constructor with particleset
   * @param target particleset
   *
   * Store mass per species and use SameMass to choose the methods.
   * if SameMass, probably faster and easy to vectorize but no impact on the performance.
   */
  BareKineticEnergy(ParticleSet& p, TrialWaveFunction& psi);
  ///destructor
  ~BareKineticEnergy() override;

  bool dependsOnWaveFunction() const override;
  std::string getClassName() const override;
  void resetTargetParticleSet(ParticleSet& p) override;

#if !defined(REMOVE_TRACEMANAGER)
  void contributeParticleQuantities() override;
  void checkoutParticleQuantities(TraceManager& tm) override;
  void deleteParticleQuantities() override;
#endif

  Return_t evaluate(ParticleSet& P) override;

  Return_t evaluateValueAndDerivatives(ParticleSet& P,
                                       const opt_variables_type& optvars,
                                       const Vector<ValueType>& dlogpsi,
                                       Vector<ValueType>& dhpsioverpsi) override;

  void mw_evaluateWithParameterDerivatives(const RefVectorWithLeader<OperatorBase>& o_list,
                                           const RefVectorWithLeader<ParticleSet>& p_list,
                                           const opt_variables_type& optvars,
                                           const RecordArray<ValueType>& dlogpsi,
                                           RecordArray<ValueType>& dhpsioverpsi) const override;

  /** Evaluate the contribution of this component for multiple walkers reporting
   *  to registered listeners from Estimators.
   */
  void mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              const std::vector<ListenerVector<RealType>>& listeners,
                              const std::vector<ListenerVector<RealType>>& ion_listeners) const override;

  /** For BareKineticEnergy since it does override any Toperator evals this needs to decay to
   *  mw_evaluatePerParticle.
   *
   *  This method must be overrideen since the default behavior is to decay to mw_evaluateWithToperator
   *  and its default behavior is to call mw_evaluate.
   */
  void mw_evaluatePerParticleWithToperator(const RefVectorWithLeader<OperatorBase>& o_list,
                                           const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                           const RefVectorWithLeader<ParticleSet>& p_list,
                                           const std::vector<ListenerVector<RealType>>& listeners,
                                           const std::vector<ListenerVector<RealType>>& ion_listeners) const override;

  /**@brief Function to compute the value, direct ionic gradient terms, and pulay terms for the local kinetic energy.
 *  
 *  This general function represents the OperatorBase interface for computing.  For an operator \hat{O}, this
 *  function will return \frac{\hat{O}\Psi_T}{\Psi_T},  \frac{\partial(\hat{O})\Psi_T}{\Psi_T}, and 
 *  \frac{\hat{O}\partial\Psi_T}{\Psi_T} - \frac{\hat{O}\Psi_T}{\Psi_T}\frac{\partial \Psi_T}{\Psi_T}.  These are 
 *  referred to as Value, HF term, and pulay term respectively.
 *
 * @param P electron particle set.
 * @param ions ion particle set
 * @param psi Trial wave function object.
 * @param hf_terms 3Nion dimensional object. All direct force terms, or ionic gradient of operator itself.
 *                 Contribution of this operator is ADDED onto hf_terms.
 * @param pulay_terms The terms coming from ionic gradients of trial wavefunction.  Contribution of this operator is
 *                 ADDED onto pulay_terms.
 * @return Value of kinetic energy operator at electron/ion positions given by P and ions.  The force contributions from
 *          this operator are added into hf_terms and pulay_terms.
 */
  Return_t evaluateWithIonDerivs(ParticleSet& P,
                                 ParticleSet& ions,
                                 TrialWaveFunction& psi,
                                 ParticleSet::ParticlePos& hf_terms,
                                 ParticleSet::ParticlePos& pulay_terms) override;

  void evaluateOneBodyOpMatrix(ParticleSet& P, const TWFFastDerivWrapper& psi, std::vector<ValueMatrix>& B) override;

  void evaluateOneBodyOpMatrixForceDeriv(ParticleSet& P,
                                         ParticleSet& source,
                                         const TWFFastDerivWrapper& psi,
                                         const int iat,
                                         std::vector<std::vector<ValueMatrix>>& Bforce) override;

#if !defined(REMOVE_TRACEMANAGER)
  Return_t evaluate_sp(ParticleSet& P);
#endif

  Return_t evaluate_orig(ParticleSet& P);

  /** implements the virtual function.
   *
   * Nothing is done but should check the mass
   */
  bool put(xmlNodePtr) override;

  bool get(std::ostream& os) const override;

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  /** initialize a shared resource and hand it to a collection
   */
  void createResource(ResourceCollection& collection) const override;

  /** acquire a shared resource from a collection
   */
  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const override;

  /** return a shared resource to a collection
   */
  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const override;

private:
  ///true, if all the species have the same mass
  bool same_mass_;

  ///mass of the particle
  FullPrecRealType particle_mass_;

  ///\f$ 1/(2 m^*) \f$
  FullPrecRealType one_over_2m_;
  ///minus_over_2m_[i] = \f$ -1/2m[i]\f$ for the ith species
  std::vector<FullPrecRealType> minus_over_2m_;

#if !defined(REMOVE_TRACEMANAGER)
  Array<TraceReal, 1>* t_sample_;
  Array<TraceComp, 1>* t_sample_comp_;
  Array<TraceComp, 2>* p_sample_;
#endif

  ParticleSet& ps_;

  struct MultiWalkerResource;
  ResourceHandle<MultiWalkerResource> mw_res_;

  TrialWaveFunction& psi_;
};

} // namespace qmcplusplus
#endif
