//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "TWFdispatcher.h"
#include "TrialWaveFunction.h"

namespace qmcplusplus
{

TWFdispatcher::TWFdispatcher(bool use_batch) : use_batch_(use_batch) {}

void TWFdispatcher::flex_evaluateLog(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                         const RefVectorWithLeader<ParticleSet>& p_list) const
{
  if (use_batch_)
    TrialWaveFunction::mw_evaluateLog(wf_list, p_list);
  else
    for(size_t iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].evaluateLog(p_list[iw]);
}

void TWFdispatcher::flex_calcRatio(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                       const RefVectorWithLeader<ParticleSet>& p_list,
                                       int iat,
                                       std::vector<PsiValueType>& ratios,
                                       ComputeType ct) const
{
  if (use_batch_)
    TrialWaveFunction::mw_calcRatio(wf_list, p_list, iat, ratios, ct);
  else
  {
    const int num_wf = wf_list.size();
    ratios.resize(num_wf);
    for(size_t iw = 0; iw < num_wf; iw++)
      ratios[iw] = wf_list[iw].calcRatio(p_list[iw], iat, ct);
  }
}

void TWFdispatcher::flex_prepareGroup(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                          const RefVectorWithLeader<ParticleSet>& p_list,
                                          int ig) const
{
  if (use_batch_)
    TrialWaveFunction::mw_prepareGroup(wf_list, p_list, ig);
  else
    for(size_t iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].prepareGroup(p_list[iw], ig);
}

void TWFdispatcher::flex_evalGrad(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                      const RefVectorWithLeader<ParticleSet>& p_list,
                                      int iat,
                                      std::vector<GradType>& grad_now) const
{
  if (use_batch_)
    TrialWaveFunction::mw_evalGrad(wf_list, p_list, iat, grad_now);
  else
  {
    const int num_wf = wf_list.size();
    grad_now.resize(num_wf);
    for(size_t iw = 0; iw < num_wf; iw++)
      grad_now[iw] = wf_list[iw].evalGrad(p_list[iw], iat);
  }
}

void TWFdispatcher::flex_calcRatioGrad(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                           const RefVectorWithLeader<ParticleSet>& p_list,
                                           int iat,
                                           std::vector<PsiValueType>& ratios,
                                           std::vector<GradType>& grad_new) const
{
  if (use_batch_)
    TrialWaveFunction::mw_calcRatioGrad(wf_list, p_list, iat, ratios, grad_new);
  else
  {
    const int num_wf = wf_list.size();
    ratios.resize(num_wf);
    grad_new.resize(num_wf);
    for(size_t iw = 0; iw < num_wf; iw++)
      ratios[iw] = wf_list[iw].calcRatioGrad(p_list[iw], iat, grad_new[iw]);
  }
}

void TWFdispatcher::flex_accept_rejectMove(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                               const RefVectorWithLeader<ParticleSet>& p_list,
                                               int iat,
                                               const std::vector<bool>& isAccepted,
                                               bool safe_to_delay) const
{
  if (use_batch_)
    TrialWaveFunction::mw_accept_rejectMove(wf_list, p_list, iat, isAccepted, safe_to_delay);
  else
    for(size_t iw = 0; iw < wf_list.size(); iw++)
      if (isAccepted[iw])
        wf_list[iw].acceptMove(p_list[iw], iat, safe_to_delay);
      else
        wf_list[iw].rejectMove(iat);
}

void TWFdispatcher::flex_completeUpdates(const RefVectorWithLeader<TrialWaveFunction>& wf_list) const
{
  if (use_batch_)
    TrialWaveFunction::mw_completeUpdates(wf_list);
  else
    for(size_t iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].completeUpdates();
}

void TWFdispatcher::flex_evaluateGL(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                        const RefVectorWithLeader<ParticleSet>& p_list,
                                        bool fromscratch) const
{
  if (use_batch_)
    TrialWaveFunction::mw_evaluateGL(wf_list, p_list, fromscratch);
  else
    for(size_t iw = 0; iw < wf_list.size(); iw++)
  {
    wf_list[iw].evaluateGL(p_list[iw], fromscratch);
    // Ye: temporal workaround to have WF.G/L always defined.
    // remove when KineticEnergy use WF.G/L instead of P.G/L
    wf_list[iw].G = p_list[iw].G;
    wf_list[iw].L = p_list[iw].L;
  }
}

void TWFdispatcher::flex_evaluateRatios(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                            const RefVector<const VirtualParticleSet>& vp_list,
                                            const RefVector<std::vector<ValueType>>& ratios_list,
                                            ComputeType ct) const
{
  if (use_batch_)
    TrialWaveFunction::mw_evaluateRatios(wf_list, vp_list, ratios_list, ct);
  else
    for(size_t iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].evaluateRatios(vp_list[iw], ratios_list[iw], ct);
}

} // namespace qmcplusplus
