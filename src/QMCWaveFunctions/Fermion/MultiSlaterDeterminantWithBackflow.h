//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MULTISLATERDETERMINANTWITHBACKFLOW_ORBITAL_H
#define QMCPLUSPLUS_MULTISLATERDETERMINANTWITHBACKFLOW_ORBITAL_H
#include <Configuration.h>
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/SPOSetProxyForMSD.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminant.h"
#include "Utilities/TimerManager.h"

namespace qmcplusplus
{
class BackflowTransformation;

/** @ingroup WaveFunctionComponent
 *  @brief MultiSlaterDeterminantWithBackflow
 */
class MultiSlaterDeterminantWithBackflow : public MultiSlaterDeterminant
{
public:
  ///constructor
  MultiSlaterDeterminantWithBackflow(ParticleSet& targetPtcl,
                                     std::vector<std::unique_ptr<SPOSetProxyForMSD>> spos,
                                     std::unique_ptr<BackflowTransformation> tr);

  ///destructor
  ~MultiSlaterDeterminantWithBackflow() override;

  void checkInVariables(opt_variables_type& active) override;
  void checkOutVariables(const opt_variables_type& active) override;
  void resetParameters(const opt_variables_type& active) override;
  void reportStatus(std::ostream& os) override;

  ValueType evaluate(const ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L) override;

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L) override;

  GradType evalGrad(ParticleSet& P, int iat) override;
  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;
  PsiValueType ratio(ParticleSet& P, int iat) override;
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;
  void restore(int iat) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;
  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override;

private:
  void resize(int, int) override;

  // transformation
  const std::unique_ptr<BackflowTransformation> BFTrans;

  // temporary storage for evaluateDerivatives
  Matrix<RealType> dpsia_up, dLa_up;
  Matrix<RealType> dpsia_dn, dLa_dn;
  Array<GradType, 3> dGa_up, dGa_dn;
};

} // namespace qmcplusplus
#endif
