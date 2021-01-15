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

void WaveFunctionComponent::mw_evaluateLog(const RefVector<WaveFunctionComponent>& WFC_list,
                                           const RefVector<ParticleSet>& P_list,
                                           const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                                           const RefVector<ParticleSet::ParticleLaplacian_t>& L_list)
{
#pragma omp parallel for
  for (int iw = 0; iw < WFC_list.size(); iw++)
    WFC_list[iw].get().evaluateLog(P_list[iw], G_list[iw], L_list[iw]);
}

void WaveFunctionComponent::mw_evalGrad(const RefVector<WaveFunctionComponent>& WFC_list,
                                        const RefVector<ParticleSet>& P_list,
                                        int iat,
                                        std::vector<GradType>& grad_now)
{
#pragma omp parallel for
  for (int iw = 0; iw < WFC_list.size(); iw++)
    grad_now[iw] = WFC_list[iw].get().evalGrad(P_list[iw].get(), iat);
}

void WaveFunctionComponent::mw_calcRatio(const RefVector<WaveFunctionComponent>& WFC_list,
                                         const RefVector<ParticleSet>& P_list,
                                         int iat,
                                         std::vector<PsiValueType>& ratios)
{
#pragma omp parallel for
  for (int iw = 0; iw < WFC_list.size(); iw++)
    ratios[iw] = WFC_list[iw].get().ratio(P_list[iw], iat);
}


PsiValueType WaveFunctionComponent::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  APP_ABORT("WaveFunctionComponent::ratioGrad is not implemented in " + ClassName + " class.");
  return ValueType();
}

void WaveFunctionComponent::ratioGradAsync(ParticleSet& P, int iat, PsiValueType& ratio, GradType& grad_iat)
{
#pragma omp task default(none) firstprivate(iat) shared(P, ratio, grad_iat)
  ratio = ratioGrad(P, iat, grad_iat);
}

void WaveFunctionComponent::mw_ratioGrad(const RefVector<WaveFunctionComponent>& WFC_list,
                                         const RefVector<ParticleSet>& P_list,
                                         int iat,
                                         std::vector<PsiValueType>& ratios,
                                         std::vector<GradType>& grad_new)
{
#pragma omp parallel for
  for (int iw = 0; iw < WFC_list.size(); iw++)
    ratios[iw] = WFC_list[iw].get().ratioGrad(P_list[iw], iat, grad_new[iw]);
}

void WaveFunctionComponent::mw_ratioGradAsync(const RefVector<WaveFunctionComponent>& WFC_list,
                                              const RefVector<ParticleSet>& P_list,
                                              int iat,
                                              std::vector<PsiValueType>& ratios,
                                              std::vector<GradType>& grad_new)
{
#pragma omp task default(none) firstprivate(WFC_list, P_list, iat) shared(ratios, grad_new)
  mw_ratioGrad(WFC_list, P_list, iat, ratios, grad_new);
}

void WaveFunctionComponent::mw_accept_rejectMove(const RefVector<WaveFunctionComponent>& WFC_list,
                                                 const RefVector<ParticleSet>& P_list,
                                                 int iat,
                                                 const std::vector<bool>& isAccepted,
                                                 bool safe_to_delay)
{
#pragma omp parallel for
  for (int iw = 0; iw < WFC_list.size(); iw++)
    if (isAccepted[iw])
      WFC_list[iw].get().acceptMove(P_list[iw], iat, safe_to_delay);
    else
      WFC_list[iw].get().restore(iat);
}

void WaveFunctionComponent::mw_completeUpdates(const RefVector<WaveFunctionComponent>& WFC_list)
{
#pragma omp parallel for
  for (int iw = 0; iw < WFC_list.size(); iw++)
    WFC_list[iw].get().completeUpdates();
}

WaveFunctionComponent::LogValueType WaveFunctionComponent::evaluateGL(ParticleSet& P,
                                                                      ParticleSet::ParticleGradient_t& G,
                                                                      ParticleSet::ParticleLaplacian_t& L,
                                                                      bool fromscratch)
{
  return evaluateLog(P, G, L);
}

void WaveFunctionComponent::mw_evaluateGL(const RefVector<WaveFunctionComponent>& WFC_list,
                                          const RefVector<ParticleSet>& P_list,
                                          const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                                          const RefVector<ParticleSet::ParticleLaplacian_t>& L_list,
                                          bool fromscratch)
{
#pragma omp parallel for
  for (int iw = 0; iw < WFC_list.size(); iw++)
    WFC_list[iw].get().evaluateGL(P_list[iw], G_list[iw], L_list[iw], fromscratch);
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

void WaveFunctionComponent::evaluateDerivRatios(VirtualParticleSet& VP,
                                                const opt_variables_type& optvars,
                                                std::vector<ValueType>& ratios,
                                                Matrix<ValueType>& dratios)
{
  //default is only ratios and zero derivatives
  evaluateRatios(VP, ratios);
}

} // namespace qmcplusplus
