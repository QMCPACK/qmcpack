//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_COUNTINGWAVEFUNCTIONCOMPONENT_H
#define QMCPLUSPLUS_COUNTINGWAVEFUNCTIONCOMPONENT_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"

#include <algorithm>
#include <memory>
#include <vector>

namespace qmcplusplus::testing
{
/** Neutral test component with clone-local observations.
 * Each clone is accessed by its owning task; inspect records after those tasks join.
 * Batched lists contain this same component type, as for other WFC batched calls.
 */
class CountingWaveFunctionComponent : public WaveFunctionComponent
{
public:
  struct ProposalRecord
  {
    bool evaluate;
    int evaluations;
  };
  struct ResolutionRecord
  {
    bool accepted;
    int evaluations;
    bool mixed;
  };

  CountingWaveFunctionComponent() = default;
  std::string getClassName() const override { return "CountingWaveFunctionComponent"; }
  LogValue evaluateLog(const ParticleSet&, ParticleSet::ParticleGradient&, ParticleSet::ParticleLaplacian&) override
  {
    ++evaluation_count_;
    log_value_ = 0;
    return log_value_;
  }
  void mw_evaluateLogAllParticles(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                  const RefVectorWithLeader<ParticleSet>& p_list,
                                  const RefVector<ParticleSet::ParticleGradient>& G_list,
                                  const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                                  const std::vector<bool>& evaluate_mask) const override
  {
    std::vector<int> before;
    before.reserve(wfc_list.size());
    for (std::size_t i = 0; i < wfc_list.size(); ++i)
      before.push_back(static_cast<const CountingWaveFunctionComponent&>(wfc_list[i]).evaluation_count_);
    WaveFunctionComponent::mw_evaluateLogAllParticles(wfc_list, p_list, G_list, L_list, evaluate_mask);
    for (std::size_t i = 0; i < wfc_list.size(); ++i)
    {
      auto& component = static_cast<CountingWaveFunctionComponent&>(wfc_list[i]);
      component.proposals_.push_back({evaluate_mask[i], component.evaluation_count_ - before[i]});
    }
  }
  void mw_accept_rejectMoveAllParticles(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                        const RefVectorWithLeader<ParticleSet>& p_list,
                                        const RefVector<ParticleSet::ParticleGradient>& G_list,
                                        const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                                        const std::vector<bool>& accepted) const override
  {
    std::vector<int> before;
    before.reserve(wfc_list.size());
    for (std::size_t i = 0; i < wfc_list.size(); ++i)
      before.push_back(static_cast<const CountingWaveFunctionComponent&>(wfc_list[i]).evaluation_count_);
    WaveFunctionComponent::mw_accept_rejectMoveAllParticles(wfc_list, p_list, G_list, L_list, accepted);
    const bool mixed = std::count(accepted.begin(), accepted.end(), true) != 0 &&
        std::count(accepted.begin(), accepted.end(), false) != 0;
    for (std::size_t i = 0; i < wfc_list.size(); ++i)
    {
      auto& component = static_cast<CountingWaveFunctionComponent&>(wfc_list[i]);
      component.resolutions_.push_back({accepted[i], component.evaluation_count_ - before[i], mixed});
    }
  }
  void acceptMove(ParticleSet&, int, bool) override {}
  void restore(int) override {}
  PsiValue ratio(ParticleSet&, int) override { return PsiValue(1); }
  void registerData(ParticleSet&, WFBufferType&) override {}
  LogValue updateBuffer(ParticleSet&, WFBufferType&, bool) override { return LogValue(0); }
  void copyFromBuffer(ParticleSet&, WFBufferType&) override {}
  void evaluateDerivatives(ParticleSet&, const OptVariables&, Vector<ValueType>&, Vector<ValueType>&) override {}
  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet&) const override
  {
    return std::make_unique<CountingWaveFunctionComponent>();
  }

  int getEvaluationCount() const { return evaluation_count_; }
  const std::vector<ProposalRecord>& getProposals() const { return proposals_; }
  const std::vector<ResolutionRecord>& getResolutions() const { return resolutions_; }

private:
  int evaluation_count_ = 0;
  std::vector<ProposalRecord> proposals_;
  std::vector<ResolutionRecord> resolutions_;
};
} // namespace qmcplusplus::testing

#endif
