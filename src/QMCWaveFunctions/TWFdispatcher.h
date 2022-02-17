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


#ifndef QMCPLUSPLUS_TWFDISPATCH_H
#define QMCPLUSPLUS_TWFDISPATCH_H

#include "TrialWaveFunction.h"
#include "TWFGrads.hpp"

namespace qmcplusplus
{
/** Wrappers for dispatching to TrialWaveFunction single walker APIs or mw_ APIs.
 * This should be only used by QMC drivers.
 * member function names must match mw_ APIs in TrialWaveFunction
 */
class TWFdispatcher
{
public:
  using PsiValueType = TrialWaveFunction::PsiValueType;
  using ComputeType  = TrialWaveFunction::ComputeType;
  using ValueType    = TrialWaveFunction::ValueType;
  using GradType     = TrialWaveFunction::GradType;
  using Complex      = TrialWaveFunction::ComplexType;

  TWFdispatcher(bool use_batch);

  void flex_evaluateLog(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                        const RefVectorWithLeader<ParticleSet>& p_list) const;

  void flex_recompute(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const std::vector<bool>& recompute) const;

  void flex_calcRatio(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      int iat,
                      std::vector<PsiValueType>& ratios,
                      ComputeType ct = ComputeType::ALL) const;

  void flex_prepareGroup(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                         const RefVectorWithLeader<ParticleSet>& p_list,
                         int ig) const;

  template<CoordsType CT>
  void flex_evalGrad(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                     const RefVectorWithLeader<ParticleSet>& p_list,
                     int iat,
                     TWFGrads<CT>& grads) const;

  template<CoordsType CT>
  void flex_calcRatioGrad(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                          const RefVectorWithLeader<ParticleSet>& p_list,
                          int iat,
                          std::vector<PsiValueType>& ratios,
                          TWFGrads<CT>& grads) const;

  void flex_accept_rejectMove(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              int iat,
                              const std::vector<bool>& isAccepted,
                              bool safe_to_delay) const;

  void flex_completeUpdates(const RefVectorWithLeader<TrialWaveFunction>& wf_list) const;

  void flex_evaluateGL(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                       const RefVectorWithLeader<ParticleSet>& p_list,
                       bool fromscratch) const;

  void flex_evaluateRatios(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                           const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                           const RefVector<std::vector<ValueType>>& ratios_list,
                           ComputeType ct) const;

private:
  bool use_batch_;
};
} // namespace qmcplusplus

#endif
