//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DIFFERENTIAL_ORBITAL_BASE_H
#define QMCPLUSPLUS_DIFFERENTIAL_ORBITAL_BASE_H
#include "Configuration.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"

/**@file DiffWaveFunctionComponent.h
 *@brief Declaration of DiffWaveFunctionComponent, NumericalDiffOrbital and AnalyticDiffOrbital
 */
namespace qmcplusplus
{
/**@defgroup WaveFunctionComponent group
 * @brief Base base class for derivatives of WaveFunctionComponent
 *
 * Derived classes implement the differentiate function which evaluates
 * - \f${\partial \log \Psi}/{\partial \alpha}\f$
 * - \f$\nabla {\partial \log \Psi}/{\partial \alpha}\f$
 * - \f$\nabla^2 {\partial \log \Psi}/{\partial \alpha}\f$
 * Each object handles one or more parameters during the optimiation.
 * The data type of refOrbital, std::vector<WaveFunctionComponent*> is intended for the cases
 * when a common variable is used by several WaveFunctionComponent class, e.g.,
 * TwoBodyJastrowOrbital<PadeFunctor> and OneBodyJastrowOrbital<PadeFunctor>
 * with one scaling parameter.
 */
struct DiffWaveFunctionComponent
{
  //@{typedefs inherited from WaveFunctionComponent
  using RealType        = WaveFunctionComponent::RealType;
  using ValueType       = WaveFunctionComponent::ValueType;
  using PosType         = WaveFunctionComponent::PosType;
  using GradVectorType  = ParticleSet::ParticleGradient;
  using ValueVectorType = ParticleSet::ParticleLaplacian;
  // the value type for \f$ log(\psi) \f$
  using LogValueType = std::complex<QMCTraits::QTFull::RealType>;
  // the value type for \f$ \psi(r')/\psi(r) \f$
  using PsiValueType = QMCTraits::QTFull::ValueType;
  //@}
  /** list of reference orbitals which contribute to the derivatives
   */
  std::vector<WaveFunctionComponent*> refOrbital;

  /// default constructor
  DiffWaveFunctionComponent(WaveFunctionComponent* orb = 0);

  ///default destructor
  virtual ~DiffWaveFunctionComponent() {}

  /** add an WaveFunctionComponent*
   */
  inline void addOrbital(WaveFunctionComponent* psi) { refOrbital.push_back(psi); }

  /** initialize the object
   *
   * Only those that need extra checks and resizing implement this function
   */
  virtual void initialize() {}

  /** evaluate derivatives at \f$\{R\}\f$
   * @param P current configuration
   * @param optvars optimizable variables
   * @param dlogpsi derivative of the log of the wavefunction
   * @param dhpsioverpsi derivative of the local kinetic energy
   */
  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   std::vector<ValueType>& dlogpsi,
                                   std::vector<ValueType>& dhpsioverpsi) = 0;

  /** evaluate derivatives at \f$\{R\}\f$
   * @param P current configuration
   * @param optvars optimizable variables
   * @param dlogpsi derivative of the log of the wavefunction
   */
  virtual void evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, std::vector<ValueType>& dlogpsi)
  {
    app_error() << "Need specialization of DiffOrbitalBase::evaluateDerivativesWF.\n";
    abort();
  }

  virtual void evaluateDerivRatios(ParticleSet& VP, const opt_variables_type& optvars, Matrix<ValueType>& dratios);

  virtual void multiplyDerivsByOrbR(std::vector<ValueType>& dlogpsi)
  {
    for (int i = 0; i < refOrbital.size(); ++i)
      for (int j = 0; j < refOrbital[i]->myVars.size(); j++)
      {
        int loc = refOrbital[j]->myVars.where(j);
        dlogpsi[loc] *= refOrbital[i]->getValue();
      }
  };

  /** check out optimizable variables
   */
  virtual void checkOutVariables(const opt_variables_type& optvars) = 0;

  /** reset the parameters during optimizations
   */
  virtual void resetParameters(const opt_variables_type& optvars) = 0;

  /** make clone
    * @param tqp target Quantum ParticleSet
    */
  virtual std::unique_ptr<DiffWaveFunctionComponent> makeClone(ParticleSet& tqp) const;
};

/** a generic DiffWaveFunctionComponent using a finite-difference method for a single optimizable parameter.
 *
 * This class handles any orbital whose analytic derivatives are not implemented nor easily available.
 */
struct NumericalDiffOrbital : public DiffWaveFunctionComponent
{
  NumericalDiffOrbital(WaveFunctionComponent* orb = 0) : DiffWaveFunctionComponent(orb) {}

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override;
  void checkOutVariables(const opt_variables_type& optvars) override;
  void resetParameters(const opt_variables_type& optvars) override;

  /// \f$\nabla \partial_{\alpha} log\Psi\f$
  GradVectorType gradLogPsi, dg_p, dg_m;
  /// \f$\nabla^2 \partial_{\alpha} log\Psi\f$
  ValueVectorType lapLogPsi, dl_p, dl_m;
};

/** a generic DiffWaveFunctionComponent using analytic derivates for a single optimizable parameter
 */
struct AnalyticDiffOrbital : public DiffWaveFunctionComponent
{
  AnalyticDiffOrbital(WaveFunctionComponent* orb = 0) : DiffWaveFunctionComponent(orb) {}

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override;
  void checkOutVariables(const opt_variables_type& optvars) override;
  void resetParameters(const opt_variables_type& optvars) override;

  /// \f$\nabla \partial_{\alpha} log\Psi\f$
  GradVectorType gradLogPsi;
  /// \f$\nabla^2 \partial_{\alpha} log\Psi\f$
  ValueVectorType lapLogPsi;
  /// get the index in the variable list
  int MyIndex;
};
} // namespace qmcplusplus
#endif
