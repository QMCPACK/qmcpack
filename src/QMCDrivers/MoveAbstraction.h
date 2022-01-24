//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MOVEABSTRACTION_H
#define QMCPLUSPLUS_MOVEABSTRACTION_H

#include "QMCDriverNew.h"
#include "Particle/PSdispatcher.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/TWFdispatcher.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDriverInput.h"
#include "type_traits/RefVectorWithLeader.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"
#include "Crowd.h"

namespace qmcplusplus
{
template<QMCDriverNew::CoordsToMove COORDS>
class MoveAbstraction
{
  using Pos     = ParticleSet::PosType;
  using Scalar  = ParticleSet::Scalar_t;
  using Real    = ParticleSet::RealType;
  using Grad    = TrialWaveFunction::GradType;
  using Complex = TrialWaveFunction::ComplexType;
  using PsiV    = TrialWaveFunction::PsiValueType;

public:
  MoveAbstraction(const PSdispatcher& ps_dispatcher,
                  const TWFdispatcher& twf_dispatcher,
                  RandomGenerator& random_gen,
                  const DriftModifierBase& drift_modifier,
                  const int num_walkers,
                  const int num_particles);

  void generateDeltas();

  void setTauForGroup(const QMCDriverInput& qmcdrv_input, const Real& invmass);

  void calcForwardMoveWithDrift(const RefVectorWithLeader<TrialWaveFunction>& twfs,
                                const RefVectorWithLeader<ParticleSet>& elecs,
                                const int iat);

  void calcForwardMove(const int iat);

  void makeMove(const RefVectorWithLeader<ParticleSet>& elecs, const int iat);

  void updateGreensFunctionWithDrift(const RefVectorWithLeader<TrialWaveFunction>& twfs,
                                     const RefVectorWithLeader<ParticleSet>& elecs,
                                     Crowd& crowd,
                                     const int iat,
                                     std::vector<PsiV>& ratios,
                                     std::vector<Real>& log_gf,
                                     std::vector<Real>& log_bg);

  void updateGreensFunction(const RefVectorWithLeader<TrialWaveFunction>& twfs,
                            const RefVectorWithLeader<ParticleSet>& elecs,
                            const int iat,
                            std::vector<PsiV>& ratios);

private:
  std::vector<Pos> drifts_;
  std::vector<Pos> walker_deltas_;
  std::vector<Grad> grads_now_, grads_new_;
  std::vector<Scalar> spindrifts_;
  std::vector<Scalar> walker_spindeltas_;
  std::vector<Complex> spingrads_now_, spingrads_new_;
  const PSdispatcher& ps_dispatcher_;
  const TWFdispatcher& twf_dispatcher_;
  RandomGenerator& random_gen_;
  const DriftModifierBase& drift_modifier_;
  Real tauovermass_;
  Real oneover2tau_;
  Real sqrttau_;
  Real spintauovermass_;
  Real oneover2spintau_;
  Real sqrtspintau_;
  const int num_walkers_;
};

template<>
MoveAbstraction<QMCDriverNew::POSITIONS>::MoveAbstraction(const PSdispatcher& ps_dispatcher,
                                                          const TWFdispatcher& twf_dispatcher,
                                                          RandomGenerator& random_gen,
                                                          const DriftModifierBase& drift_modifier,
                                                          const int num_walkers,
                                                          const int num_particles)
    : ps_dispatcher_(ps_dispatcher),
      twf_dispatcher_(twf_dispatcher),
      random_gen_(random_gen),
      drift_modifier_(drift_modifier),
      num_walkers_(num_walkers)
{
  drifts_.resize(num_walkers_);
  walker_deltas_.resize(num_walkers_ * num_particles);
  grads_now_.resize(num_walkers_);
  grads_new_.resize(num_walkers_);
}

template<>
MoveAbstraction<QMCDriverNew::POSITIONS_SPINS>::MoveAbstraction(const PSdispatcher& ps_dispatcher,
                                                                const TWFdispatcher& twf_dispatcher,
                                                                RandomGenerator& random_gen,
                                                                const DriftModifierBase& drift_modifier,
                                                                const int num_walkers,
                                                                const int num_particles)
    : ps_dispatcher_(ps_dispatcher),
      twf_dispatcher_(twf_dispatcher),
      random_gen_(random_gen),
      drift_modifier_(drift_modifier),
      num_walkers_(num_walkers)
{
  drifts_.resize(num_walkers_);
  walker_deltas_.resize(num_walkers_ * num_particles);
  grads_now_.resize(num_walkers_);
  grads_new_.resize(num_walkers_);
  spindrifts_.resize(num_walkers_);
  walker_spindeltas_.resize(num_walkers_ * num_particles);
  spingrads_now_.resize(num_walkers_);
  spingrads_new_.resize(num_walkers_);
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS>::generateDeltas()
{
  makeGaussRandomWithEngine(walker_deltas_, random_gen_);
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS_SPINS>::generateDeltas()
{
  makeGaussRandomWithEngine(walker_deltas_, random_gen_);
  makeGaussRandomWithEngine(walker_spindeltas_, random_gen_);
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS>::setTauForGroup(const QMCDriverInput& qmcdrv_input, const Real& invmass)
{
  tauovermass_ = qmcdrv_input.get_tau() * invmass;
  oneover2tau_ = 0.5 / tauovermass_;
  sqrttau_     = std::sqrt(tauovermass_);
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS_SPINS>::setTauForGroup(const QMCDriverInput& qmcdrv_input,
                                                                    const Real& invmass)
{
  tauovermass_     = qmcdrv_input.get_tau() * invmass;
  oneover2tau_     = 0.5 / tauovermass_;
  sqrttau_         = std::sqrt(tauovermass_);
  spintauovermass_ = tauovermass_ / qmcdrv_input.get_spin_mass();
  oneover2spintau_ = 0.5 / spintauovermass_;
  sqrtspintau_     = std::sqrt(spintauovermass_);
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS>::calcForwardMoveWithDrift(
    const RefVectorWithLeader<TrialWaveFunction>& twfs,
    const RefVectorWithLeader<ParticleSet>& elecs,
    const int iat)
{
  auto delta_r_start = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;

  twf_dispatcher_.flex_evalGrad(twfs, elecs, iat, grads_now_);
  drift_modifier_.getDrifts(tauovermass_, grads_now_, drifts_);
  std::transform(drifts_.begin(), drifts_.end(), delta_r_start, drifts_.begin(),
                 [st = sqrttau_](const Pos& drift, const Pos& delta_r) { return drift + (st * delta_r); });
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS_SPINS>::calcForwardMoveWithDrift(
    const RefVectorWithLeader<TrialWaveFunction>& twfs,
    const RefVectorWithLeader<ParticleSet>& elecs,
    const int iat)
{
  auto delta_r_start    = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end      = delta_r_start + num_walkers_;
  auto delta_spin_start = walker_spindeltas_.begin() + iat * num_walkers_;
  auto delta_spin_end   = delta_spin_start + num_walkers_;

  twf_dispatcher_.flex_evalGradWithSpin(twfs, elecs, iat, grads_now_, spingrads_now_);
  drift_modifier_.getDrifts(tauovermass_, grads_now_, drifts_);
  std::transform(drifts_.begin(), drifts_.end(), delta_r_start, drifts_.begin(),
                 [st = sqrttau_](const Pos& drift, const Pos& delta_r) { return drift + (st * delta_r); });
  //spin part of forward move
  drift_modifier_.getDrifts(spintauovermass_, spingrads_now_, spindrifts_);
  std::transform(spindrifts_.begin(), spindrifts_.end(), delta_spin_start, spindrifts_.begin(),
                 [st = sqrtspintau_](const Scalar& spindrift, const Scalar& delta_spin) {
                   return spindrift + (st * delta_spin);
                 });
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS>::calcForwardMove(const int iat)
{
  auto delta_r_start = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;

  std::transform(delta_r_start, delta_r_end, drifts_.begin(),
                 [st = sqrttau_](const Pos& delta_r) { return st * delta_r; });
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS_SPINS>::calcForwardMove(const int iat)
{
  auto delta_r_start    = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end      = delta_r_start + num_walkers_;
  auto delta_spin_start = walker_spindeltas_.begin() + iat * num_walkers_;
  auto delta_spin_end   = delta_spin_start + num_walkers_;

  std::transform(delta_r_start, delta_r_end, drifts_.begin(),
                 [st = sqrttau_](const Pos& delta_r) { return st * delta_r; });
  std::transform(delta_spin_start, delta_spin_end, spindrifts_.begin(),
                 [st = sqrtspintau_](const Scalar& delta_spin) { return st * delta_spin; });
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS>::makeMove(const RefVectorWithLeader<ParticleSet>& elecs, const int iat)
{
  ps_dispatcher_.flex_makeMove(elecs, iat, drifts_);
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS_SPINS>::makeMove(const RefVectorWithLeader<ParticleSet>& elecs,
                                                              const int iat)
{
  ps_dispatcher_.flex_makeMoveWithSpin(elecs, iat, drifts_, spindrifts_);
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS>::updateGreensFunctionWithDrift(
    const RefVectorWithLeader<TrialWaveFunction>& twfs,
    const RefVectorWithLeader<ParticleSet>& elecs,
    Crowd& crowd,
    const int iat,
    std::vector<PsiV>& ratios,
    std::vector<Real>& log_gf,
    std::vector<Real>& log_gb)
{
  auto delta_r_start = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;

  twf_dispatcher_.flex_calcRatioGrad(twfs, elecs, iat, ratios, grads_new_);
  std::transform(delta_r_start, delta_r_end, log_gf.begin(), [](const Pos& delta_r) {
    constexpr Real mhalf(-0.5);
    return mhalf * dot(delta_r, delta_r);
  });
  drift_modifier_.getDrifts(tauovermass_, grads_new_, drifts_);

  std::transform(crowd.beginElectrons(), crowd.endElectrons(), drifts_.begin(), drifts_.begin(),
                 [iat](const ParticleSet& ps, const Pos& drift) { return ps.R[iat] - ps.getActivePos() - drift; });

  std::transform(drifts_.begin(), drifts_.end(), log_gb.begin(),
                 [halfovertau = oneover2tau_](const Pos& drift) { return -halfovertau * dot(drift, drift); });
}

template<>
void MoveAbstraction<QMCDriverNew::POSITIONS_SPINS>::updateGreensFunctionWithDrift(
    const RefVectorWithLeader<TrialWaveFunction>& twfs,
    const RefVectorWithLeader<ParticleSet>& elecs,
    Crowd& crowd,
    const int iat,
    std::vector<PsiV>& ratios,
    std::vector<Real>& log_gf,
    std::vector<Real>& log_gb)
{
  auto delta_r_start    = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end      = delta_r_start + num_walkers_;
  auto delta_spin_start = walker_spindeltas_.begin() + iat * num_walkers_;
  auto delta_spin_end   = delta_spin_start + num_walkers_;

  twf_dispatcher_.flex_calcRatioGradWithSpin(twfs, elecs, iat, ratios, grads_new_, spingrads_new_);
  std::transform(delta_r_start, delta_r_end, log_gf.begin(), [](const Pos& delta_r) {
    constexpr Real mhalf(-0.5);
    return mhalf * dot(delta_r, delta_r);
  });
  //add spin part to greens function
  std::transform(delta_spin_start, delta_spin_end, log_gf.begin(), log_gf.begin(),
                 [](const Scalar& delta_spin, const Real& loggf) {
                   constexpr Real mhalf(-0.5);
                   return loggf + mhalf * delta_spin * delta_spin;
                 });

  drift_modifier_.getDrifts(tauovermass_, grads_new_, drifts_);
  drift_modifier_.getDrifts(spintauovermass_, spingrads_new_, spindrifts_);

  std::transform(crowd.beginElectrons(), crowd.endElectrons(), drifts_.begin(), drifts_.begin(),
                 [iat](const ParticleSet& ps, const Pos& drift) { return ps.R[iat] - ps.getActivePos() - drift; });
  std::transform(crowd.beginElectrons(), crowd.endElectrons(), spindrifts_.begin(), spindrifts_.begin(),
                 [iat](const ParticleSet& ps, const Scalar& spindrift) {
                   return ps.spins[iat] - ps.getActiveSpinVal() - spindrift;
                 });

  std::transform(drifts_.begin(), drifts_.end(), log_gb.begin(),
                 [halfovertau = oneover2tau_](const Pos& drift) { return -halfovertau * dot(drift, drift); });
  //add spin part to greens function
  std::transform(spindrifts_.begin(), spindrifts_.end(), log_gb.begin(), log_gb.begin(),
                 [halfovertau = oneover2spintau_](const Scalar& spindrift, const Real& loggb) {
                   return loggb - halfovertau * spindrift * spindrift;
                 });
}
template<QMCDriverNew::CoordsToMove COORDS>
void MoveAbstraction<COORDS>::updateGreensFunction(const RefVectorWithLeader<TrialWaveFunction>& twfs,
                                                   const RefVectorWithLeader<ParticleSet>& elecs,
                                                   const int iat,
                                                   std::vector<PsiV>& ratios)
{
  twf_dispatcher_.flex_calcRatio(twfs, elecs, iat, ratios);
}

} // namespace qmcplusplus

#endif
