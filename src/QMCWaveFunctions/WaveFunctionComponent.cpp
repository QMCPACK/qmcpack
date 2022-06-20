//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "WaveFunctionComponent.h"

namespace qmcplusplus
{
// for return types
using PsiValueType = WaveFunctionComponent::PsiValueType;

WaveFunctionComponent::WaveFunctionComponent(const std::string& class_name, const std::string& obj_name)
    : IsOptimizing(false),
      Optimizable(true),
      is_fermionic(false),
      UpdateMode(ORB_WALKER),
      ClassName(class_name),
      myName(obj_name),
      Bytes_in_WFBuffer(0),
      log_value_(0.0)
{
  if (ClassName.empty())
    throw std::runtime_error("WaveFunctionComponent ClassName cannot be empty!");
}

WaveFunctionComponent::~WaveFunctionComponent() = default;

void WaveFunctionComponent::mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                           const RefVectorWithLeader<ParticleSet>& p_list,
                                           const RefVector<ParticleSet::ParticleGradient>& G_list,
                                           const RefVector<ParticleSet::ParticleLaplacian>& L_list) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list[iw].evaluateLog(p_list[iw], G_list[iw], L_list[iw]);
}

void WaveFunctionComponent::recompute(const ParticleSet& P)
{
  ParticleSet::ParticleGradient temp_G(P.getTotalNum());
  ParticleSet::ParticleLaplacian temp_L(P.getTotalNum());

  evaluateLog(P, temp_G, temp_L);
}

void WaveFunctionComponent::mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                         const RefVectorWithLeader<ParticleSet>& p_list,
                                         const std::vector<bool>& recompute) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    if (recompute[iw])
      wfc_list[iw].recompute(p_list[iw]);
}

void WaveFunctionComponent::mw_prepareGroup(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                            const RefVectorWithLeader<ParticleSet>& p_list,
                                            int ig) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list[iw].prepareGroup(p_list[iw], ig);
}

template<CoordsType CT>
void WaveFunctionComponent::mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                        const RefVectorWithLeader<ParticleSet>& p_list,
                                        const int iat,
                                        TWFGrads<CT>& grad_now) const
{
  if constexpr (CT == CoordsType::POS_SPIN)
    mw_evalGradWithSpin(wfc_list, p_list, iat, grad_now.grads_positions, grad_now.grads_spins);
  else
    mw_evalGrad(wfc_list, p_list, iat, grad_now.grads_positions);
}

void WaveFunctionComponent::mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                        const RefVectorWithLeader<ParticleSet>& p_list,
                                        int iat,
                                        std::vector<GradType>& grad_now) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    grad_now[iw] = wfc_list[iw].evalGrad(p_list[iw], iat);
}

void WaveFunctionComponent::mw_evalGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                const RefVectorWithLeader<ParticleSet>& p_list,
                                                int iat,
                                                std::vector<GradType>& grad_now,
                                                std::vector<ComplexType>& spingrad_now) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    grad_now[iw] = wfc_list[iw].evalGradWithSpin(p_list[iw], iat, spingrad_now[iw]);
}

void WaveFunctionComponent::mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                         const RefVectorWithLeader<ParticleSet>& p_list,
                                         int iat,
                                         std::vector<PsiValueType>& ratios) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    ratios[iw] = wfc_list[iw].ratio(p_list[iw], iat);
}


PsiValueType WaveFunctionComponent::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  APP_ABORT("WaveFunctionComponent::ratioGrad is not implemented in " + ClassName + " class.");
  return ValueType();
}

template<CoordsType CT>
void WaveFunctionComponent::mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                         const RefVectorWithLeader<ParticleSet>& p_list,
                                         int iat,
                                         std::vector<PsiValueType>& ratios,
                                         TWFGrads<CT>& grad_new) const
{
  if constexpr (CT == CoordsType::POS_SPIN)
    mw_ratioGradWithSpin(wfc_list, p_list, iat, ratios, grad_new.grads_positions, grad_new.grads_spins);
  else
    mw_ratioGrad(wfc_list, p_list, iat, ratios, grad_new.grads_positions);
}

void WaveFunctionComponent::mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                         const RefVectorWithLeader<ParticleSet>& p_list,
                                         int iat,
                                         std::vector<PsiValueType>& ratios,
                                         std::vector<GradType>& grad_new) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    ratios[iw] = wfc_list[iw].ratioGrad(p_list[iw], iat, grad_new[iw]);
}

void WaveFunctionComponent::mw_ratioGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                 const RefVectorWithLeader<ParticleSet>& p_list,
                                                 int iat,
                                                 std::vector<PsiValueType>& ratios,
                                                 std::vector<GradType>& grad_new,
                                                 std::vector<ComplexType>& spingrad_new) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    ratios[iw] = wfc_list[iw].ratioGradWithSpin(p_list[iw], iat, grad_new[iw], spingrad_new[iw]);
}

void WaveFunctionComponent::mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                 const RefVectorWithLeader<ParticleSet>& p_list,
                                                 int iat,
                                                 const std::vector<bool>& isAccepted,
                                                 bool safe_to_delay) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    if (isAccepted[iw])
      wfc_list[iw].acceptMove(p_list[iw], iat, safe_to_delay);
    else
      wfc_list[iw].restore(iat);
}

void WaveFunctionComponent::mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list[iw].completeUpdates();
}

WaveFunctionComponent::LogValueType WaveFunctionComponent::evaluateGL(const ParticleSet& P,
                                                                      ParticleSet::ParticleGradient& G,
                                                                      ParticleSet::ParticleLaplacian& L,
                                                                      bool fromscratch)
{
  return evaluateLog(P, G, L);
}

void WaveFunctionComponent::mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                          const RefVectorWithLeader<ParticleSet>& p_list,
                                          const RefVector<ParticleSet::ParticleGradient>& G_list,
                                          const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                                          bool fromscratch) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list[iw].evaluateGL(p_list[iw], G_list[iw], L_list[iw], fromscratch);
}

void WaveFunctionComponent::evaluateDerivativesWF(ParticleSet& P,
                                                  const opt_variables_type& active,
                                                  std::vector<ValueType>& dlogpsi)
{
  throw std::runtime_error("WaveFunctionComponent::evaluateDerivativesWF is not implemented by " + ClassName);
}

/*@todo makeClone should be a pure virtual function
 */
std::unique_ptr<WaveFunctionComponent> WaveFunctionComponent::makeClone(ParticleSet& tpq) const
{
  APP_ABORT("Implement WaveFunctionComponent::makeClone " + ClassName + " class.");
  return std::unique_ptr<WaveFunctionComponent>();
}

WaveFunctionComponent::RealType WaveFunctionComponent::KECorrection() { return 0; }

void WaveFunctionComponent::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  assert(P.getTotalNum() == ratios.size());
  for (int i = 0; i < P.getTotalNum(); ++i)
    ratios[i] = ratio(P, i);
}

void WaveFunctionComponent::evaluateRatios(const VirtualParticleSet& P, std::vector<ValueType>& ratios)
{
  std::ostringstream o;
  o << "WaveFunctionComponent::evaluateRatios is not implemented by " << ClassName;
  APP_ABORT(o.str());
}

void WaveFunctionComponent::mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                              const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                              std::vector<std::vector<ValueType>>& ratios) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list[iw].evaluateRatios(vp_list[iw], ratios[iw]);
}

void WaveFunctionComponent::evaluateDerivRatios(const VirtualParticleSet& VP,
                                                const opt_variables_type& optvars,
                                                std::vector<ValueType>& ratios,
                                                Matrix<ValueType>& dratios)
{
  //default is only ratios and zero derivatives
  evaluateRatios(VP, ratios);
}

void WaveFunctionComponent::registerTWFFastDerivWrapper(const ParticleSet& P, TWFFastDerivWrapper& twf) const
{
  std::ostringstream o;
  o << "WaveFunctionComponent::registerTWFFastDerivWrapper is not implemented by " << ClassName;
  APP_ABORT(o.str());
}

template void WaveFunctionComponent::mw_evalGrad<CoordsType::POS>(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int iat,
    TWFGrads<CoordsType::POS>& grad_now) const;
template void WaveFunctionComponent::mw_evalGrad<CoordsType::POS_SPIN>(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int iat,
    TWFGrads<CoordsType::POS_SPIN>& grad_now) const;
template void WaveFunctionComponent::mw_ratioGrad<CoordsType::POS>(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int iat,
    std::vector<PsiValueType>& ratios,
    TWFGrads<CoordsType::POS>& grad_new) const;
template void WaveFunctionComponent::mw_ratioGrad<CoordsType::POS_SPIN>(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int iat,
    std::vector<PsiValueType>& ratios,
    TWFGrads<CoordsType::POS_SPIN>& grad_new) const;

} // namespace qmcplusplus
