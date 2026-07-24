//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/DeepQMC/DeepQMCWaveFunctionComponent.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>

namespace qmcplusplus
{
namespace
{
void validateResultShape(const DeepQMCBridge::BatchResult& result, int batch_size, int n_elec)
{
  const std::size_t batch     = static_cast<std::size_t>(batch_size);
  const std::size_t particles = static_cast<std::size_t>(batch_size) * static_cast<std::size_t>(n_elec);
  if (result.log_values.size() != batch)
    throw std::runtime_error("DeepQMC bridge returned wrong number of log values");
  if (result.grad_log_values.size() != particles * OHMMS_DIM)
    throw std::runtime_error("DeepQMC bridge returned wrong number of gradient values");
  if (result.lap_log_values.size() != particles)
    throw std::runtime_error("DeepQMC bridge returned wrong number of laplacian values");
}

} // namespace

DeepQMCWaveFunctionComponent::DeepQMCWaveFunctionComponent(std::string name,
                                                           const ParticleSet& ions,
                                                           std::shared_ptr<const DeepQMCBridge> bridge,
                                                           int mol_idx)
    : WaveFunctionComponent(std::move(name)), ions_(ions), bridge_(std::move(bridge)), mol_idx_(mol_idx)
{
  if (!bridge_)
    throw std::runtime_error("DeepQMCWaveFunctionComponent requires a non-null DeepQMCBridge");
}

std::vector<DeepQMCWaveFunctionComponent::RealType> DeepQMCWaveFunctionComponent::flattenIonCoords(
    const ParticleSet& ions)
{
  std::vector<RealType> ion_coords;
  ion_coords.reserve(ions.getTotalNum() * OHMMS_DIM);
  for (int iat = 0; iat < ions.getTotalNum(); ++iat)
    for (int d = 0; d < OHMMS_DIM; ++d)
      ion_coords.push_back(ions.R[iat][d]);
  return ion_coords;
}

void DeepQMCWaveFunctionComponent::appendElectronCoords(const ParticleSet& electrons,
                                                        std::vector<RealType>& electron_coords,
                                                        int active_iat)
{
  for (int iat = 0; iat < electrons.getTotalNum(); ++iat)
  {
    const auto& pos = (iat == active_iat) ? electrons.activeR(iat) : electrons.R[iat];
    for (int d = 0; d < OHMMS_DIM; ++d)
      electron_coords.push_back(pos[d]);
  }
}

DeepQMCBridge::BatchResult DeepQMCWaveFunctionComponent::evaluateBatch(
    const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    int active_iat)
{
  if (wfc_list.size() != p_list.size())
    throw std::runtime_error("DeepQMCWaveFunctionComponent::evaluateBatch list size mismatch");
  if (wfc_list.empty())
    return {};

  const auto& leader   = wfc_list.getCastedLeader<DeepQMCWaveFunctionComponent>();
  const int batch_size = static_cast<int>(wfc_list.size());
  const int n_elec     = static_cast<int>(p_list[0].getTotalNum());

  std::vector<RealType> electron_coords;
  electron_coords.reserve(static_cast<std::size_t>(batch_size) * static_cast<std::size_t>(n_elec) * OHMMS_DIM);
  for (int iw = 0; iw < batch_size; ++iw)
  {
    if (p_list[iw].getTotalNum() != n_elec)
      throw std::runtime_error(
          "DeepQMCWaveFunctionComponent requires all walkers in a batch to have the same electron count");
    appendElectronCoords(p_list[iw], electron_coords, active_iat);
  }

  const std::vector<RealType> ion_coords = flattenIonCoords(leader.ions_);
  DeepQMCBridge::BatchResult result =
      leader.bridge_->evaluateLogBatch(ion_coords, electron_coords, leader.mol_idx_, batch_size, n_elec);
  validateResultShape(result, batch_size, n_elec);
  return result;
}

DeepQMCBridge::BatchResult DeepQMCWaveFunctionComponent::evaluateOne(const ParticleSet& electrons, int active_iat) const
{
  auto& mutable_p = const_cast<ParticleSet&>(electrons);
  RefVectorWithLeader<WaveFunctionComponent> wfc_list(*const_cast<DeepQMCWaveFunctionComponent*>(this));
  wfc_list.push_back(*const_cast<DeepQMCWaveFunctionComponent*>(this));
  RefVectorWithLeader<ParticleSet> p_list(mutable_p);
  p_list.push_back(mutable_p);
  return evaluateBatch(wfc_list, p_list, active_iat);
}

DeepQMCWaveFunctionComponent::LogValue DeepQMCWaveFunctionComponent::evaluateLog(const ParticleSet& P,
                                                                                 ParticleSet::ParticleGradient& G,
                                                                                 ParticleSet::ParticleLaplacian& L)
{
  auto& mutable_p = const_cast<ParticleSet&>(P);
  RefVectorWithLeader<WaveFunctionComponent> wfc_list(*this);
  wfc_list.push_back(*this);
  RefVectorWithLeader<ParticleSet> p_list(mutable_p);
  p_list.push_back(mutable_p);
  RefVector<ParticleSet::ParticleGradient> G_list;
  G_list.push_back(G);
  RefVector<ParticleSet::ParticleLaplacian> L_list;
  L_list.push_back(L);

  mw_evaluateLog(wfc_list, p_list, G_list, L_list);
  return log_value_;
}

void DeepQMCWaveFunctionComponent::mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                  const RefVectorWithLeader<ParticleSet>& p_list,
                                                  const RefVector<ParticleSet::ParticleGradient>& G_list,
                                                  const RefVector<ParticleSet::ParticleLaplacian>& L_list) const
{
  if (wfc_list.size() != p_list.size() || wfc_list.size() != G_list.size() || wfc_list.size() != L_list.size())
    throw std::runtime_error("DeepQMCWaveFunctionComponent::mw_evaluateLog list size mismatch");
  if (wfc_list.empty())
    return;

  const int batch_size              = static_cast<int>(wfc_list.size());
  const int n_elec                  = static_cast<int>(p_list[0].getTotalNum());
  DeepQMCBridge::BatchResult result = evaluateBatch(wfc_list, p_list);

  for (int iw = 0; iw < batch_size; ++iw)
  {
    auto& component      = wfc_list.getCastedElement<DeepQMCWaveFunctionComponent>(iw);
    component.log_value_ = LogValue(result.log_values[iw]);

    auto& G = G_list[iw].get();
    auto& L = L_list[iw].get();
    if (G.size() < n_elec)
      G.resize(n_elec);
    if (L.size() < n_elec)
      L.resize(n_elec);

    for (int iat = 0; iat < n_elec; ++iat)
    {
      const std::size_t particle_offset = (static_cast<std::size_t>(iw) * n_elec + iat);
      GradType grad;
      for (int d = 0; d < OHMMS_DIM; ++d)
        grad[d] = result.grad_log_values[particle_offset * OHMMS_DIM + d];
      G[iat] += grad;
      L[iat] += result.lap_log_values[particle_offset];
    }
  }
}

void DeepQMCWaveFunctionComponent::mw_prepareGroup(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                                   int ig) const
{}

void DeepQMCWaveFunctionComponent::mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                               const RefVectorWithLeader<ParticleSet>& p_list,
                                               int iat,
                                               std::vector<GradType>& grad_now) const
{
  if (wfc_list.size() != p_list.size() || wfc_list.size() != grad_now.size())
    throw std::runtime_error("DeepQMCWaveFunctionComponent::mw_evalGrad list size mismatch");
  if (wfc_list.empty())
    return;

  const int n_elec  = p_list[0].getTotalNum();
  const auto result = evaluateBatch(wfc_list, p_list);
  for (int iw = 0; iw < wfc_list.size(); ++iw)
    for (int d = 0; d < OHMMS_DIM; ++d)
      grad_now[iw][d] = result.grad_log_values[(static_cast<std::size_t>(iw) * n_elec + iat) * OHMMS_DIM + d];
}

void DeepQMCWaveFunctionComponent::mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                const RefVectorWithLeader<ParticleSet>& p_list,
                                                int iat,
                                                std::vector<PsiValue>& ratios) const
{
  if (wfc_list.size() != p_list.size())
    throw std::runtime_error("DeepQMCWaveFunctionComponent::mw_calcRatio list size mismatch");
  const int batch_size = static_cast<int>(wfc_list.size());
  ratios.resize(batch_size);
  if (wfc_list.empty())
    return;

  const auto result = evaluateBatch(wfc_list, p_list, iat);
  for (int iw = 0; iw < batch_size; ++iw)
  {
    auto& component                   = wfc_list.getCastedElement<DeepQMCWaveFunctionComponent>(iw);
    component.proposed_log_value_     = LogValue(result.log_values[iw]);
    component.has_proposed_log_value_ = true;
    ratios[iw]                        = std::exp(std::real(component.proposed_log_value_ - component.log_value_));
  }
}

void DeepQMCWaveFunctionComponent::mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                const RefVectorWithLeader<ParticleSet>& p_list,
                                                int iat,
                                                std::vector<PsiValue>& ratios,
                                                std::vector<GradType>& grad_new) const
{
  if (wfc_list.size() != p_list.size() || wfc_list.size() != grad_new.size())
    throw std::runtime_error("DeepQMCWaveFunctionComponent::mw_ratioGrad list size mismatch");
  const int batch_size = static_cast<int>(wfc_list.size());
  ratios.resize(batch_size);
  if (wfc_list.empty())
    return;

  const int n_elec  = p_list[0].getTotalNum();
  const auto result = evaluateBatch(wfc_list, p_list, iat);
  for (int iw = 0; iw < batch_size; ++iw)
  {
    auto& component                   = wfc_list.getCastedElement<DeepQMCWaveFunctionComponent>(iw);
    component.proposed_log_value_     = LogValue(result.log_values[iw]);
    component.has_proposed_log_value_ = true;
    ratios[iw]                        = std::exp(std::real(component.proposed_log_value_ - component.log_value_));
    for (int d = 0; d < OHMMS_DIM; ++d)
      grad_new[iw][d] = result.grad_log_values[(static_cast<std::size_t>(iw) * n_elec + iat) * OHMMS_DIM + d];
  }
}

void DeepQMCWaveFunctionComponent::mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                        const RefVectorWithLeader<ParticleSet>& p_list,
                                                        int iat,
                                                        const std::vector<bool>& isAccepted,
                                                        bool safe_to_delay) const
{
  if (wfc_list.size() != isAccepted.size())
    throw std::runtime_error("DeepQMCWaveFunctionComponent::mw_accept_rejectMove list size mismatch");

  for (int iw = 0; iw < wfc_list.size(); ++iw)
  {
    auto& component = wfc_list.getCastedElement<DeepQMCWaveFunctionComponent>(iw);
    if (isAccepted[iw] && component.has_proposed_log_value_)
      component.log_value_ = component.proposed_log_value_;
    component.has_proposed_log_value_ = false;
  }
}

void DeepQMCWaveFunctionComponent::mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{}

void DeepQMCWaveFunctionComponent::mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                 const RefVectorWithLeader<ParticleSet>& p_list,
                                                 const RefVector<ParticleSet::ParticleGradient>& G_list,
                                                 const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                                                 bool fromscratch) const
{
  mw_evaluateLog(wfc_list, p_list, G_list, L_list);
}

void DeepQMCWaveFunctionComponent::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  if (has_proposed_log_value_)
    log_value_ = proposed_log_value_;
  has_proposed_log_value_ = false;
}

void DeepQMCWaveFunctionComponent::restore(int iat) { has_proposed_log_value_ = false; }

DeepQMCWaveFunctionComponent::PsiValue DeepQMCWaveFunctionComponent::ratio(ParticleSet& P, int iat)
{
  const auto result       = evaluateOne(P, iat);
  proposed_log_value_     = LogValue(result.log_values[0]);
  has_proposed_log_value_ = true;
  return std::exp(std::real(proposed_log_value_ - log_value_));
}

DeepQMCWaveFunctionComponent::GradType DeepQMCWaveFunctionComponent::evalGrad(ParticleSet& P, int iat)
{
  const auto result = evaluateOne(P);
  GradType grad;
  for (int d = 0; d < OHMMS_DIM; ++d)
    grad[d] = result.grad_log_values[static_cast<std::size_t>(iat) * OHMMS_DIM + d];
  return grad;
}

DeepQMCWaveFunctionComponent::PsiValue DeepQMCWaveFunctionComponent::ratioGrad(ParticleSet& P,
                                                                               int iat,
                                                                               GradType& grad_iat)
{
  const auto result       = evaluateOne(P, iat);
  proposed_log_value_     = LogValue(result.log_values[0]);
  has_proposed_log_value_ = true;
  for (int d = 0; d < OHMMS_DIM; ++d)
    grad_iat[d] = result.grad_log_values[static_cast<std::size_t>(iat) * OHMMS_DIM + d];
  return std::exp(std::real(proposed_log_value_ - log_value_));
}

DeepQMCWaveFunctionComponent::LogValue DeepQMCWaveFunctionComponent::updateBuffer(ParticleSet& P,
                                                                                  WFBufferType& buf,
                                                                                  bool fromscratch)
{
  return evaluateLog(P, P.G, P.L);
}

void DeepQMCWaveFunctionComponent::evaluateDerivatives(ParticleSet& P,
                                                       const OptVariables& optvars,
                                                       Vector<ValueType>& dlogpsi,
                                                       Vector<ValueType>& dhpsioverpsi)
{
  // Pretrained DeepQMC inference parameters are fixed in this prototype.
}

std::unique_ptr<WaveFunctionComponent> DeepQMCWaveFunctionComponent::makeClone(ParticleSet& tpq) const
{
  return std::make_unique<DeepQMCWaveFunctionComponent>(my_name_, ions_, bridge_, mol_idx_);
}

} // namespace qmcplusplus
