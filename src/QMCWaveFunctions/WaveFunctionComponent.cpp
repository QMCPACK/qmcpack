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
#include "QMCWaveFunctions/DiffWaveFunctionComponent.h"

namespace qmcplusplus
{
// for return types
using PsiValueType = WaveFunctionComponent::PsiValueType;

WaveFunctionComponent::WaveFunctionComponent(const std::string& class_name, const std::string& obj_name)
    : IsOptimizing(false),
      Optimizable(true),
      is_fermionic(false),
      UpdateMode(ORB_WALKER),
      LogValue(0.0),
      dPsi(nullptr),
      ClassName(class_name),
      myName(obj_name),
      Bytes_in_WFBuffer(0)
{
  if (ClassName.empty())
    throw std::runtime_error("WaveFunctionComponent ClassName cannot be empty!");
}

void WaveFunctionComponent::mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                           const RefVectorWithLeader<ParticleSet>& p_list,
                                           const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                                           const RefVector<ParticleSet::ParticleLaplacian_t>& L_list) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list[iw].evaluateLog(p_list[iw], G_list[iw], L_list[iw]);
}

void WaveFunctionComponent::recompute(const ParticleSet& P)
{
  ParticleSet::ParticleGradient_t temp_G(P.getTotalNum());
  ParticleSet::ParticleLaplacian_t temp_L(P.getTotalNum());

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
                                                                      ParticleSet::ParticleGradient_t& G,
                                                                      ParticleSet::ParticleLaplacian_t& L,
                                                                      bool fromscratch)
{
  return evaluateLog(P, G, L);
}

void WaveFunctionComponent::mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                          const RefVectorWithLeader<ParticleSet>& p_list,
                                          const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                                          const RefVector<ParticleSet::ParticleLaplacian_t>& L_list,
                                          bool fromscratch) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list[iw].evaluateGL(p_list[iw], G_list[iw], L_list[iw], fromscratch);
}

void WaveFunctionComponent::setDiffOrbital(DiffWaveFunctionComponentPtr d) { dPsi = d; }

void WaveFunctionComponent::evaluateDerivatives(ParticleSet& P,
                                                const opt_variables_type& active,
                                                std::vector<ValueType>& dlogpsi,
                                                std::vector<ValueType>& dhpsioverpsi)
{
  if (dPsi)
    dPsi->evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi);
}

void WaveFunctionComponent::evaluateDerivativesWF(ParticleSet& P,
                                                  const opt_variables_type& active,
                                                  std::vector<ValueType>& dlogpsi)
{
  if (dPsi)
    dPsi->evaluateDerivativesWF(P, active, dlogpsi);
}

/*@todo makeClone should be a pure virtual function
 */
WaveFunctionComponentPtr WaveFunctionComponent::makeClone(ParticleSet& tpq) const
{
  APP_ABORT("Implement WaveFunctionComponent::makeClone " + ClassName + " class.");
  return 0;
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
                                              const RefVector<const VirtualParticleSet>& vp_list,
                                              std::vector<std::vector<ValueType>>& ratios) const
{
  assert(this == &wfc_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list[iw].evaluateRatios(vp_list[iw], ratios[iw]);
}

void WaveFunctionComponent::evaluateDerivRatios(VirtualParticleSet& VP,
                                                const opt_variables_type& optvars,
                                                std::vector<ValueType>& ratios,
                                                Matrix<ValueType>& dratios)
{
  //default is only ratios and zero derivatives
  evaluateRatios(VP, ratios);
}

} // namespace qmcplusplus
