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
  MCCoords<CT> drifts_;
  /// all of the spatial diffusion moves for all walkers and particles for a given step
  MCCoords<CT> walker_deltas_;
  /// spatial gradients for a single electron across multiple walkers
  std::vector<Grad> grads_now_, grads_new_;
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
  std::unique_ptr<Taus<Real, CT>> taus_;
  const int num_walkers_;
};

template<CoordsType CT>
inline MoveAbstraction<CT>::MoveAbstraction(const PSdispatcher& ps_dispatcher,
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
  if constexpr (CT == CoordsType::POS_SPIN)
  {
    spingrads_now_.resize(num_walkers_);
    spingrads_new_.resize(num_walkers_);
  }
}

template<CoordsType CT>
inline void MoveAbstraction<CT>::generateDeltas()
{
  makeGaussRandomWithEngine(walker_deltas_.positions, random_gen_);
  if constexpr (CT == CoordsType::POS_SPIN)
    makeGaussRandomWithEngine(walker_deltas_.spins, random_gen_);
}

template<CoordsType CT>
inline void MoveAbstraction<CT>::setTauForGroup(const QMCDriverInput& qmcdrv_input, const Real& invmass)
{
  if constexpr (CT == CoordsType::POS_SPIN)
    taus_ = std::make_unique<Taus<Real, CoordsType::POS_SPIN>>(qmcdrv_input.get_tau(), invmass,
                                                               qmcdrv_input.get_spin_mass());
  else
    taus_ = std::make_unique<Taus<Real, CT>>(qmcdrv_input.get_tau(), invmass);
}

template<CoordsType CT>
inline void MoveAbstraction<CT>::calcForwardMoveWithDrift(const TWFdispatcher& twf_dispatcher,
                                                          const RefVectorWithLeader<TrialWaveFunction>& twfs,
                                                          const int iat)
{
  auto delta_r_start = walker_deltas_.positions.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;

  twf_dispatcher.flex_evalGrad(twfs, elecs_, iat, grads_now_);
  drift_modifier_.getDrifts(taus_->tauovermass, grads_now_, drifts_.positions);
  std::transform(drifts_.positions.begin(), drifts_.positions.end(), delta_r_start, drifts_.positions.begin(),
                 [st = taus_->sqrttau](const Pos& drift, const Pos& delta_r) { return drift + (st * delta_r); });

  if constexpr (CT == CoordsType::POS_SPIN)
  {
    auto delta_spin_start = walker_deltas_.spins.begin() + iat * num_walkers_;
    auto delta_spin_end   = delta_spin_start + num_walkers_;
    drift_modifier_.getDrifts(taus_->spin_tauovermass, spingrads_now_, drifts_.spins);
    std::transform(drifts_.spins.begin(), drifts_.spins.end(), delta_spin_start, drifts_.spins.begin(),
                   [st = taus_->spin_sqrttau](const Scalar& spindrift, const Scalar& delta_spin) {
                     return spindrift + (st * delta_spin);
                   });
  }
}

template<CoordsType CT>
inline void MoveAbstraction<CT>::calcForwardMove(const int iat)
{
  auto delta_r_start = walker_deltas_.positions.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;

  std::transform(delta_r_start, delta_r_end, drifts_.positions.begin(),
                 [st = taus_->sqrttau](const Pos& delta_r) { return st * delta_r; });

  if constexpr (CT == CoordsType::POS_SPIN)
  {
    auto delta_spin_start = walker_deltas_.spins.begin() + iat * num_walkers_;
    auto delta_spin_end   = delta_spin_start + num_walkers_;
    std::transform(delta_spin_start, delta_spin_end, drifts_.spins.begin(),
                   [st = taus_->spin_sqrttau](const Scalar& delta_spin) { return st * delta_spin; });
  }
}

template<CoordsType CT>
inline void MoveAbstraction<CT>::makeMove(const int iat)
{
  ps_dispatcher_.flex_makeMove(elecs_, iat, drifts_);
}


template<CoordsType CT>
inline void MoveAbstraction<CT>::updateGreensFunctionWithDrift(const TWFdispatcher& twf_dispatcher,
                                                               const RefVectorWithLeader<TrialWaveFunction>& twfs,
                                                               const int iat,
                                                               std::vector<PsiV>& ratios,
                                                               std::vector<Real>& log_gf,
                                                               std::vector<Real>& log_gb)
{
  if constexpr (CT == CoordsType::POS_SPIN)
    twf_dispatcher.flex_calcRatioGradWithSpin(twfs, elecs_, iat, ratios, grads_new_, spingrads_new_);
  else
    twf_dispatcher.flex_calcRatioGrad(twfs, elecs_, iat, ratios, grads_new_);

  auto delta_r_start = walker_deltas_.positions.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;

  std::transform(delta_r_start, delta_r_end, log_gf.begin(), [](const Pos& delta_r) {
    constexpr Real mhalf(-0.5);
    return mhalf * dot(delta_r, delta_r);
  });

  drift_modifier_.getDrifts(taus_->tauovermass, grads_new_, drifts_.positions);

  std::transform(elecs_.begin(), elecs_.end(), drifts_.positions.begin(), drifts_.positions.begin(),
                 [iat](const ParticleSet& ps, const Pos& drift) { return ps.R[iat] - ps.getActivePos() - drift; });

  std::transform(drifts_.positions.begin(), drifts_.positions.end(), log_gb.begin(),
                 [halfovertau = taus_->oneover2tau](const Pos& drift) { return -halfovertau * dot(drift, drift); });

  if constexpr (CT == CoordsType::POS_SPIN)
  {
    auto delta_spin_start = walker_deltas_.spins.begin() + iat * num_walkers_;
    auto delta_spin_end   = delta_spin_start + num_walkers_;
    std::transform(delta_spin_start, delta_spin_end, log_gf.begin(), log_gf.begin(),
                   [](const Scalar& delta_spin, const Real& loggf) {
                     constexpr Real mhalf(-0.5);
                     return loggf + mhalf * delta_spin * delta_spin;
                   });
    drift_modifier_.getDrifts(taus_->spin_tauovermass, spingrads_new_, drifts_.spins);
    std::transform(elecs_.begin(), elecs_.end(), drifts_.spins.begin(), drifts_.spins.begin(),
                   [iat](const ParticleSet& ps, const Scalar& spindrift) {
                     return ps.spins[iat] - ps.getActiveSpinVal() - spindrift;
                   });
    std::transform(drifts_.spins.begin(), drifts_.spins.end(), log_gb.begin(), log_gb.begin(),
                   [halfovertau = taus_->spin_oneover2tau](const Scalar& spindrift, const Real& loggb) {
                     return loggb - halfovertau * spindrift * spindrift;
                   });
  }
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
  auto delta_r_start = walker_deltas_.positions.begin() + iat * num_walkers_;
  auto delta_r_end   = delta_r_start + num_walkers_;
  assert(rr.size() == delta_r_end - delta_r_start);
  std::transform(delta_r_start, delta_r_end, rr.begin(),
                 [t = taus_->tauovermass](auto& delta_r) { return t * dot(delta_r, delta_r); });
}

} // namespace qmcplusplus

#endif
