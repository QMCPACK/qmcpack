//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "TWFdispatcher.h"
#include <cassert>
#include "TrialWaveFunction.h"

namespace qmcplusplus
{
TWFdispatcher::TWFdispatcher(bool use_batch) : use_batch_(use_batch) {}

void TWFdispatcher::flex_evaluateLog(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                     const RefVectorWithLeader<ParticleSet>& p_list) const
{
  assert(wf_list.size() == p_list.size());
  if (use_batch_)
    TrialWaveFunction::mw_evaluateLog(wf_list, p_list);
  else
    for (size_t iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].evaluateLog(p_list[iw]);
}

void TWFdispatcher::flex_recompute(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                   const std::vector<bool>& recompute) const
{
  assert(wf_list.size() == p_list.size());
  if (use_batch_)
    TrialWaveFunction::mw_recompute(wf_list, p_list, recompute);
  else
    for (size_t iw = 0; iw < wf_list.size(); iw++)
      if (recompute[iw])
        wf_list[iw].recompute(p_list[iw]);
}

void TWFdispatcher::flex_calcRatio(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                   int iat,
                                   std::vector<PsiValueType>& ratios,
                                   ComputeType ct) const
{
  assert(wf_list.size() == p_list.size());
  if (use_batch_)
    TrialWaveFunction::mw_calcRatio(wf_list, p_list, iat, ratios, ct);
  else
  {
    const int num_wf = wf_list.size();
    ratios.resize(num_wf);
    for (size_t iw = 0; iw < num_wf; iw++)
      ratios[iw] = wf_list[iw].calcRatio(p_list[iw], iat, ct);
  }
}

void TWFdispatcher::flex_prepareGroup(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                      const RefVectorWithLeader<ParticleSet>& p_list,
                                      int ig) const
{
  assert(wf_list.size() == p_list.size());
  if (use_batch_)
    TrialWaveFunction::mw_prepareGroup(wf_list, p_list, ig);
  else
    for (size_t iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].prepareGroup(p_list[iw], ig);
}

template<CoordsType CT>
void TWFdispatcher::flex_evalGrad(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                  const RefVectorWithLeader<ParticleSet>& p_list,
                                  int iat,
                                  TWFGrads<CT>& grads) const
{
  assert(wf_list.size() == p_list.size());
  if (use_batch_)
    TrialWaveFunction::mw_evalGrad(wf_list, p_list, iat, grads);
  else
  {
    const int num_wf = wf_list.size();
    grads.resize(num_wf);
    for (size_t iw = 0; iw < num_wf; iw++)
      if constexpr (CT == CoordsType::POS_SPIN)
        grads.grads_positions[iw] = wf_list[iw].evalGradWithSpin(p_list[iw], iat, grads.grads_spins[iw]);
      else
        grads.grads_positions[iw] = wf_list[iw].evalGrad(p_list[iw], iat);
  }
}

template<CoordsType CT>
void TWFdispatcher::flex_calcRatioGrad(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                       const RefVectorWithLeader<ParticleSet>& p_list,
                                       int iat,
                                       std::vector<PsiValueType>& ratios,
                                       TWFGrads<CT>& grads) const
{
  assert(wf_list.size() == p_list.size());
  if (use_batch_)
    TrialWaveFunction::mw_calcRatioGrad(wf_list, p_list, iat, ratios, grads);
  else
  {
    const int num_wf = wf_list.size();
    ratios.resize(num_wf);
    grads.resize(num_wf);
    for (size_t iw = 0; iw < num_wf; iw++)
      if constexpr (CT == CoordsType::POS_SPIN)
        ratios[iw] = wf_list[iw].calcRatioGradWithSpin(p_list[iw], iat, grads.grads_positions[iw], grads.grads_spins[iw]);
      else
        ratios[iw] = wf_list[iw].calcRatioGrad(p_list[iw], iat, grads.grads_positions[iw]);
  }
}

void TWFdispatcher::flex_accept_rejectMove(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                           const RefVectorWithLeader<ParticleSet>& p_list,
                                           int iat,
                                           const std::vector<bool>& isAccepted,
                                           bool safe_to_delay) const
{
  assert(wf_list.size() == p_list.size());
  if (use_batch_)
    TrialWaveFunction::mw_accept_rejectMove(wf_list, p_list, iat, isAccepted, safe_to_delay);
  else
    for (size_t iw = 0; iw < wf_list.size(); iw++)
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
    for (TrialWaveFunction& wf : wf_list)
      wf.completeUpdates();
}

void TWFdispatcher::flex_evaluateGL(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                    const RefVectorWithLeader<ParticleSet>& p_list,
                                    bool fromscratch) const
{
  assert(wf_list.size() == p_list.size());
  if (use_batch_)
    TrialWaveFunction::mw_evaluateGL(wf_list, p_list, fromscratch);
  else
    for (size_t iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].evaluateGL(p_list[iw], fromscratch);
}

void TWFdispatcher::flex_evaluateRatios(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                        const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                        const RefVector<std::vector<ValueType>>& ratios_list,
                                        ComputeType ct) const
{
  assert(wf_list.size() == vp_list.size());
  assert(wf_list.size() == ratios_list.size());
  if (use_batch_)
    TrialWaveFunction::mw_evaluateRatios(wf_list, vp_list, ratios_list, ct);
  else
    for (size_t iw = 0; iw < wf_list.size(); iw++)
      wf_list[iw].evaluateRatios(vp_list[iw], ratios_list[iw], ct);
}

template void TWFdispatcher::flex_evalGrad<CoordsType::POS>(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                                            const RefVectorWithLeader<ParticleSet>& p_list,
                                                            int iat,
                                                            TWFGrads<CoordsType::POS>& grads) const;
template void TWFdispatcher::flex_evalGrad<CoordsType::POS_SPIN>(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                                                 const RefVectorWithLeader<ParticleSet>& p_list,
                                                                 int iat,
                                                                 TWFGrads<CoordsType::POS_SPIN>& grads) const;
template void TWFdispatcher::flex_calcRatioGrad<CoordsType::POS>(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                                                 const RefVectorWithLeader<ParticleSet>& p_list,
                                                                 int iat,
                                                                 std::vector<PsiValueType>& ratios,
                                                                 TWFGrads<CoordsType::POS>& grads) const;
template void TWFdispatcher::flex_calcRatioGrad<CoordsType::POS_SPIN>(
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int iat,
    std::vector<PsiValueType>& ratios,
    TWFGrads<CoordsType::POS_SPIN>& grads) const;

} // namespace qmcplusplus
