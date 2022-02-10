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
  ///true, if all the species have the same mass
  bool SameMass;
  ///mass of the particle
  FullPrecRealType M;
  ///\f$ 1/(2 m^*) \f$
  FullPrecRealType OneOver2M;
  ///MinusOver2M[i] = \f$ -1/2m[i]\f$ for the ith species
  std::vector<FullPrecRealType> MinusOver2M;

  ParticleSet::ParticleGradient Gtmp;
  ParticleSet::ParticleLaplacian Ltmp;

  ///single particle trace samples
  bool streaming_kinetic;
  bool streaming_kinetic_comp;
  bool streaming_momentum;

#if !defined(REMOVE_TRACEMANAGER)
  Array<TraceReal, 1>* T_sample;
  Array<TraceComp, 1>* T_sample_comp;
  Array<TraceComp, 2>* p_sample;
#endif
  ParticleSet& Ps;

  /** constructor with particleset
   * @param target particleset
   *
   * Store mass per species and use SameMass to choose the methods.
   * if SameMass, probably faster and easy to vectorize but no impact on the performance.
   */
  BareKineticEnergy(ParticleSet& p);
  ///destructor
  ~BareKineticEnergy() override;

  void resetTargetParticleSet(ParticleSet& P) override {}

#if !defined(REMOVE_TRACEMANAGER)
  void contributeParticleQuantities() override;
  void checkoutParticleQuantities(TraceManager& tm) override;
  void deleteParticleQuantities() override;
#endif

  Return_t evaluate(ParticleSet& P) override;

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
                                         const ParticleSet& source,
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

#ifdef QMC_CUDA
  ////////////////////////////////
  // Vectorized version for GPU //
  ////////////////////////////////
  // Nothing is done on GPU here, just copy into vector
  void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy) override;
#endif
};

} // namespace qmcplusplus
#endif
