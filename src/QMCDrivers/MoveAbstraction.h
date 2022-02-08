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

#include "Particle/PSdispatcher.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/TWFdispatcher.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDriverInput.h"
#include "type_traits/RefVectorWithLeader.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"
#include "Particle/MCCoords.hpp"

namespace qmcplusplus
{
/** abstraction class to handle particle moves in batched QMC drivers
 *
 * Templated on CoordsType defined in MC
 * Currently supports CoordsType::POS and CoordsType::POS_SPIN, which includes dynamic spins in particle moves
 */
template<CoordsType CT>
class MoveAbstraction
{
  using Pos     = ParticleSet::PosType;
  using Scalar  = ParticleSet::Scalar_t;
  using Real    = ParticleSet::RealType;
  using Grad    = TrialWaveFunction::GradType;
  using Complex = TrialWaveFunction::ComplexType;
  using PsiV    = TrialWaveFunction::PsiValueType;

public:
  /** Constructor
   *
   * uses references to dispatchers, rng, drift_modifers, etc which are used regularly for the particle moves
   * resizes all of the vectors to be used throughout this move
   * \param[in] random_gen
   */
  MoveAbstraction(const PSdispatcher& ps_dispatcher,
                  const RefVectorWithLeader<ParticleSet>& elecs,
                  RandomGenerator& random_gen,
                  const DriftModifierBase& drift_modifier,
                  const int num_walkers,
                  const int num_particles);

  /** generates an entire steps worth of displacements
   *
   * Uses gaussian random numbers to generate the displacements for all coordinates over all walkers
   */
  void generateDeltas();

  /** sets the timestep information for the move
   *
   * e.g. spatial only moves (CoordsType::POS) only need the spatial timestep,  whereas
   * spin moves need a spin timestep defined by the spin mass provided by the driver input
   */
  void setTauForGroup(const QMCDriverInput& qmcdrv_input, const Real& invmass);

  /** Calulates the forward move for all particle coordinates
   * 
   * updates all the walkers in the crowd for particle iat using 
   * \f[
   * \mathbf{r}' - \mathbf{r} = \sqrt{\tau} \eta + \tau \mathbf{v}
   * \f] 
   * using drift and diffusion. 
   * This also includes the spin drift and diffusion depending on the template parameter
   */
  void calcForwardMoveWithDrift(const TWFdispatcher& twf_dispatcher,
                                const RefVectorWithLeader<TrialWaveFunction>& twfs,
                                const int iat);

  /** Calculates the forward move for all particle coordinates without drift
   *
   * updates all the walkers in the crowd for particle iat using 
   * \f[
   * \mathbf{r}' - \mathbf{r} = \sqrt{\tau}\eta 
   * \f]
   * i.e., diffusion only. Only used in VMCBatched
   * This also includes the spin diffusion depending on the template parameter.
   */
  void calcForwardMove(const int iat);

  /** makes particle move
   *
   * updates particle iat for all walkers in the crowd. This uses the PSdispatcher to use flex_ move APIs
   */
  void makeMove(const int iat);

  /** calculates the greens functions for forward and reverse moves and TWF ratio, used for acceptance in QMCDrivers
   *
   * includes the drift part of the move, i.e.
   * \f[
   * G(\mathbf{r}'\leftarrow\mathbf{r}, \tau = \exp\left[-\frac{1}{2\tau}|\mathbf{r}'-\mathbf{r} - \tau \mathbf{v}(\mathbf{r})|^2 \right]
   * \f]
   * and the ratio
   * \f[
   * \frac{\Psi_T(\mathbf{r}')}{\Psi_T(\mathbf{r})}
   * \f]
   * This is the necessary data to calculate the acceptance ratio in the QMCDriver
   * also adds to the spin move component depending on the template parameter
   * \param[out] ratios ratio of trial wave functions for all walkers for particle iat
   * \param[out] log_gf log of greens function for forward move for all walkers for particle iat
   * \param[out] log_gb log of greens function for reverse move for all walkers for particle iat
   */
  void updateGreensFunctionWithDrift(const TWFdispatcher& twf_dispatcher,
                                     const RefVectorWithLeader<TrialWaveFunction>& twfs,
                                     const int iat,
                                     std::vector<PsiV>& ratios,
                                     std::vector<Real>& log_gf,
                                     std::vector<Real>& log_bg);

  /** calculates TWF ratio, used for acceptance in QMCDrivers
   *
   * does not include the drift part of the move.
   * \f[
   * G(\mathbf{r}'\leftarrow\mathbf{r}, \tau = \exp\left[-\frac{1}{2\tau}|\mathbf{r}'-\mathbf{r}|^2 \right]
   * \f]
   * Therefore, in the acceptance ratio this cancels since \f$G(\mathbf{r}'\leftarrow \mathbf{r}, \tau) = G(\mathbf{r}\leftarrow\mathbf{r}',\tau)$\f. 
   * Therefore, we only need to calculate the ratios
   * \f[
   * \frac{\Psi_T(\mathbf{r}')}{\Psi_T(\mathbf{r})}
   * \f]
   * This is the necessary data to calculate the acceptance ratio in the VMCBatched, without drift
   * \param[out] ratios ratio of trial wave functions for all walkers for particle iat
   */
  void updateGreensFunction(const TWFdispatcher& twf_dispatcher,
                            const RefVectorWithLeader<TrialWaveFunction>& twfs,
                            const int iat,
                            std::vector<PsiV>& ratios);

  /** accumulates the data to  help construct the effective timestep
   *
   * the effective timestep used in DMC algorithm is 
   * \f[
   * \tau_{\rm eff} = \frac{\sum R_{\rm accepted}^2}{\sum R_{\rm proposed}^2}
   * \f]
   * which is accumulated over all electrons in a walker
   *
   * rr below is the \f$r^2$\f for the current particle iat for the walkers and will be accumulated for
   * each particle in the DMC driver.
   * \param[out] rr
   */
  void updaterr(const int iat, std::vector<Real>& rr);

private:
  /// spatial drift part of move for a single particle across multiple walkers
  std::vector<Pos> drifts_;
  /// all of the spatial diffusion moves for all walkers and particles for a given step
  std::vector<Pos> walker_deltas_;
  /// spatial gradients for a single electron across multiple walkers
  std::vector<Grad> grads_now_, grads_new_;
  /// spin drift part of move for a single particle across multiple walkers
  std::vector<Scalar> spindrifts_;
  /// all of the spin diffusion moves for all walkerss and particles for a given step
  std::vector<Scalar> walker_spindeltas_;
  /// spin gradients for a single electrons across multiple walkers
  std::vector<Complex> spingrads_now_, spingrads_new_;
  /// provides flex_ APIs to do select sw/mw updates of particle set
  const PSdispatcher& ps_dispatcher_;
  /// ParticleSets for each walker
  const RefVectorWithLeader<ParticleSet>& elecs_;
  /// rng, provided by ContextForSteps
  RandomGenerator& random_gen_;
  /// drift modifer used to limit size of drift velocity
  const DriftModifierBase& drift_modifier_;
  /// timesteps
  std::unique_ptr<Taus<Real,CT>> taus_;
  const int num_walkers_;
};

template<>
inline MoveAbstraction<CoordsType::POS>::MoveAbstraction(const PSdispatcher& ps_dispatcher,
                                                         const RefVectorWithLeader<ParticleSet>& elecs,
                                                         RandomGenerator& random_gen,
                                                         const DriftModifierBase& drift_modifier,
                                                         const int num_walkers,
                                                         const int num_particles)
    : ps_dispatcher_(ps_dispatcher),
      elecs_(elecs),
      random_gen_(random_gen),
      drift_modifier_(drift_modifier),
      taus_(nullptr),
      num_walkers_(num_walkers)
{
  drifts_.resize(num_walkers_);
  walker_deltas_.resize(num_walkers_ * num_particles);
  grads_now_.resize(num_walkers_);
  grads_new_.resize(num_walkers_);
}

template<>
inline MoveAbstraction<CoordsType::POS_SPIN>::MoveAbstraction(const PSdispatcher& ps_dispatcher,
                                                              const RefVectorWithLeader<ParticleSet>& elecs,
                                                              RandomGenerator& random_gen,
                                                              const DriftModifierBase& drift_modifier,
                                                              const int num_walkers,
                                                              const int num_particles)
    : ps_dispatcher_(ps_dispatcher),
      elecs_(elecs),
      random_gen_(random_gen),
      drift_modifier_(drift_modifier),
      taus_(nullptr),
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
inline void MoveAbstraction<CoordsType::POS>::generateDeltas()
{
  makeGaussRandomWithEngine(walker_deltas_, random_gen_);
}

template<>
inline void MoveAbstraction<CoordsType::POS_SPIN>::generateDeltas()
{
  makeGaussRandomWithEngine(walker_deltas_, random_gen_);
  makeGaussRandomWithEngine(walker_spindeltas_, random_gen_);
}

template<>
inline void MoveAbstraction<CoordsType::POS>::setTauForGroup(const QMCDriverInput& qmcdrv_input, const Real& invmass)
{
    taus_ = std::make_unique<Taus<Real,CoordsType::POS>>(qmcdrv_input.get_tau(), invmass);
}

template<>
inline void MoveAbstraction<CoordsType::POS_SPIN>::setTauForGroup(const QMCDriverInput& qmcdrv_input,
                                                                  const Real& invmass)
{
    taus_ = std::make_unique<Taus<Real,CoordsType::POS_SPIN>>(qmcdrv_input.get_tau(), invmass, qmcdrv_input.get_spin_mass());
}

template<>
inline void MoveAbstraction<CoordsType::POS>::calcForwardMoveWithDrift(
    const TWFdispatcher& twf_dispatcher,
    const RefVectorWithLeader<TrialWaveFunction>& twfs,
    const int iat)
{
  auto delta_r_start = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;

  twf_dispatcher.flex_evalGrad(twfs, elecs_, iat, grads_now_);
  drift_modifier_.getDrifts(taus_->tauovermass, grads_now_, drifts_);
  std::transform(drifts_.begin(), drifts_.end(), delta_r_start, drifts_.begin(),
                 [st = taus_->sqrttau](const Pos& drift, const Pos& delta_r) { return drift + (st * delta_r); });
}

template<>
inline void MoveAbstraction<CoordsType::POS_SPIN>::calcForwardMoveWithDrift(
    const TWFdispatcher& twf_dispatcher,
    const RefVectorWithLeader<TrialWaveFunction>& twfs,
    const int iat)
{
  auto delta_r_start    = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end      = delta_r_start + num_walkers_;
  auto delta_spin_start = walker_spindeltas_.begin() + iat * num_walkers_;
  auto delta_spin_end   = delta_spin_start + num_walkers_;

  twf_dispatcher.flex_evalGradWithSpin(twfs, elecs_, iat, grads_now_, spingrads_now_);
  drift_modifier_.getDrifts(taus_->tauovermass, grads_now_, drifts_);
  std::transform(drifts_.begin(), drifts_.end(), delta_r_start, drifts_.begin(),
                 [st = taus_->sqrttau](const Pos& drift, const Pos& delta_r) { return drift + (st * delta_r); });
  //spin part of forward move
  drift_modifier_.getDrifts(taus_->spin_tauovermass, spingrads_now_, spindrifts_);
  std::transform(spindrifts_.begin(), spindrifts_.end(), delta_spin_start, spindrifts_.begin(),
                 [st = taus_->spin_sqrttau](const Scalar& spindrift, const Scalar& delta_spin) {
                   return spindrift + (st * delta_spin);
                 });
}

template<>
inline void MoveAbstraction<CoordsType::POS>::calcForwardMove(const int iat)
{
  auto delta_r_start = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;

  std::transform(delta_r_start, delta_r_end, drifts_.begin(),
                 [st = taus_->sqrttau](const Pos& delta_r) { return st * delta_r; });
}

template<>
inline void MoveAbstraction<CoordsType::POS_SPIN>::calcForwardMove(const int iat)
{
  auto delta_r_start    = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end      = delta_r_start + num_walkers_;
  auto delta_spin_start = walker_spindeltas_.begin() + iat * num_walkers_;
  auto delta_spin_end   = delta_spin_start + num_walkers_;

  std::transform(delta_r_start, delta_r_end, drifts_.begin(),
                 [st = taus_->sqrttau](const Pos& delta_r) { return st * delta_r; });
  std::transform(delta_spin_start, delta_spin_end, spindrifts_.begin(),
                 [st = taus_->spin_sqrttau](const Scalar& delta_spin) { return st * delta_spin; });
}

template<>
inline void MoveAbstraction<CoordsType::POS>::makeMove(const int iat)
{
  ps_dispatcher_.flex_makeMove(elecs_, iat, drifts_);
}

template<>
inline void MoveAbstraction<CoordsType::POS_SPIN>::makeMove(const int iat)
{
  ps_dispatcher_.flex_makeMove(elecs_, iat, drifts_);
  ParticleSet& elec_leader = elecs_.getLeader();
  elec_leader.mw_makeSpinMove(elecs_, iat, spindrifts_);
}

template<>
inline void MoveAbstraction<CoordsType::POS>::updateGreensFunctionWithDrift(
    const TWFdispatcher& twf_dispatcher,
    const RefVectorWithLeader<TrialWaveFunction>& twfs,
    const int iat,
    std::vector<PsiV>& ratios,
    std::vector<Real>& log_gf,
    std::vector<Real>& log_gb)
{
  auto delta_r_start = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;

  twf_dispatcher.flex_calcRatioGrad(twfs, elecs_, iat, ratios, grads_new_);
  std::transform(delta_r_start, delta_r_end, log_gf.begin(), [](const Pos& delta_r) {
    constexpr Real mhalf(-0.5);
    return mhalf * dot(delta_r, delta_r);
  });
  drift_modifier_.getDrifts(taus_->tauovermass, grads_new_, drifts_);

  std::transform(elecs_.begin(), elecs_.end(), drifts_.begin(), drifts_.begin(),
                 [iat](const ParticleSet& ps, const Pos& drift) { return ps.R[iat] - ps.getActivePos() - drift; });

  std::transform(drifts_.begin(), drifts_.end(), log_gb.begin(),
                 [halfovertau = taus_->oneover2tau](const Pos& drift) { return -halfovertau * dot(drift, drift); });
}

template<>
inline void MoveAbstraction<CoordsType::POS_SPIN>::updateGreensFunctionWithDrift(
    const TWFdispatcher& twf_dispatcher,
    const RefVectorWithLeader<TrialWaveFunction>& twfs,
    const int iat,
    std::vector<PsiV>& ratios,
    std::vector<Real>& log_gf,
    std::vector<Real>& log_gb)
{
  auto delta_r_start    = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end      = delta_r_start + num_walkers_;
  auto delta_spin_start = walker_spindeltas_.begin() + iat * num_walkers_;
  auto delta_spin_end   = delta_spin_start + num_walkers_;

  twf_dispatcher.flex_calcRatioGradWithSpin(twfs, elecs_, iat, ratios, grads_new_, spingrads_new_);
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

  drift_modifier_.getDrifts(taus_->tauovermass, grads_new_, drifts_);
  drift_modifier_.getDrifts(taus_->spin_tauovermass, spingrads_new_, spindrifts_);

  std::transform(elecs_.begin(), elecs_.end(), drifts_.begin(), drifts_.begin(),
                 [iat](const ParticleSet& ps, const Pos& drift) { return ps.R[iat] - ps.getActivePos() - drift; });
  std::transform(elecs_.begin(), elecs_.end(), spindrifts_.begin(), spindrifts_.begin(),
                 [iat](const ParticleSet& ps, const Scalar& spindrift) {
                   return ps.spins[iat] - ps.getActiveSpinVal() - spindrift;
                 });

  std::transform(drifts_.begin(), drifts_.end(), log_gb.begin(),
                 [halfovertau = taus_->oneover2tau](const Pos& drift) { return -halfovertau * dot(drift, drift); });
  //add spin part to greens function
  std::transform(spindrifts_.begin(), spindrifts_.end(), log_gb.begin(), log_gb.begin(),
                 [halfovertau = taus_->spin_oneover2tau](const Scalar& spindrift, const Real& loggb) {
                   return loggb - halfovertau * spindrift * spindrift;
                 });
}

template<CoordsType CT>
inline void MoveAbstraction<CT>::updateGreensFunction(const TWFdispatcher& twf_dispatcher,
                                                      const RefVectorWithLeader<TrialWaveFunction>& twfs,
                                                      const int iat,
                                                      std::vector<PsiV>& ratios)
{
  twf_dispatcher.flex_calcRatio(twfs, elecs_, iat, ratios);
}

template<CoordsType CT>
inline void MoveAbstraction<CT>::updaterr(const int iat, std::vector<Real>& rr)
{
  auto delta_r_start = walker_deltas_.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;
  assert(rr.size() == delta_r_end - delta_r_start);
  std::transform(delta_r_start, delta_r_end, rr.begin(),
                 [t = taus_->tauovermass](auto& delta_r) { return t * dot(delta_r, delta_r); });
}

} // namespace qmcplusplus

#endif
