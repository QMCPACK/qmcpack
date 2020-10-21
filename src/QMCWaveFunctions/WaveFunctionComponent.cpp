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

// WaveFunctionComponent::WaveFunctionComponent(const WaveFunctionComponent& old):
//   Optimizable(old.Optimizable), UseBuffer(old.UseBuffer),
//   dPsi(old.dPsi),dLogPsi(old.dLogPsi),d2LogPsi(old.d2LogPsi),
//   ClassName(old.ClassName),myVars(old.myVars)
// {
//   //
//   //if(dLogPsi.size()) dLogPsi.resize(dLogPsi.size());
//   //if(d2LogPsi.size()) dLogPsi.resize(d2LogPsi.size());
//   //if(dPsi) dPsi=old.dPsi->makeClone();
// }


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
