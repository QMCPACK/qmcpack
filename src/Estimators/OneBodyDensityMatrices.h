//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Based on code from: QMCHamiltonians/DensityMatrices1B.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ONE_BODY_DENSITY_MATRICES_H
#define QMCPLUSPLUS_ONE_BODY_DENSITY_MATRICES_H

#include <vector>
#include <functional>

#include "Estimators/OperatorEstBase.h"
#include "type_traits/complex_help.hpp"
#include "QMCWaveFunctions/CompositeSPOSet.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "OneBodyDensityMatricesInput.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include <SpeciesSet.h>
#include <StdRandom.h>

namespace qmcplusplus
{

namespace testing
{
template<typename T>
class OneBodyDensityMatricesTests;
}

/** Per crowd Estimator for OneBodyDensityMatrices aka 1RDM DensityMatrices1B
 *
 *  \todo most matrices are written to by incrementing a single vector index
 *        into their memory. This isn't compatible with aligned rows and ignores
 *        much less error prone accessing that Matrix supplys.  Fix it.
 *  \todo functions favor output arguments or state updates over return values,
 *        simplify.
 */
class OneBodyDensityMatrices : public OperatorEstBase
{
public:
  using Value         = QMCTraits::ValueType;
  using FullPrecValue = QMCTraits::FullPrecValueType;
  using Real          = RealAlias<Value>;
  using FullPrecReal  = RealAlias<FullPrecValue>;
  using Grad          = TinyVector<Value, OHMMS_DIM>;
  using Lattice       = PtclOnLatticeTraits::ParticleLayout;
  using Position      = QMCTraits::PosType;

  using Evaluator  = OneBodyDensityMatricesInput::Evaluator;
  using Integrator = OneBodyDensityMatricesInput::Integrator;

  enum class Sampling
  {
    VOLUME_BASED,
    METROPOLIS,
    NO_SAMPLING
  };

private:
  OneBodyDensityMatricesInput input_;
  Lattice lattice_;
  SpeciesSet species_;

  /** @ingroup Derived simulation parameters determined by computation based in input
   *  @{
   */
  /// samples_ are altered based on integrator_ so always have a simulation object copy.
  size_t samples_;
  /// Sampling method, this derived values in input_
  Sampling sampling_;
  /// If not defined in OBDMI taken from lattice_
  Position center_;
  /// with respect to center_ using lattice_;
  Position rcorner_;
  Real volume_;
  bool periodic_;
  /** @} */

  //data members \todo analyze lifecycles allocation optimization or state?
  CompositeSPOSet basis_functions_;
  Vector<Value> basis_values_;
  Vector<Value> basis_norms_;
  Vector<Grad> basis_gradients_;
  Vector<Value> basis_laplacians_;
  std::vector<Position> rsamples_;
  Vector<Real> samples_weights_;
  int basis_size_;
  std::vector<int> species_sizes_;
  std::vector<std::string> species_names_;
  /** @ingroup Working space, I'm assuming not longterm state.
   *  @{ */
  /** per particle ratios
   *  size: particles
   */
  std::vector<Value> psi_ratios_;

  /// row major per sample workspaces
  /** conj(basis_values) for each particle 
   *  size: samples * basis_size
   *  vector is over species
   *  each matrix row: particle column: basis_value
   */
  std::vector<Matrix<Value>> Phi_NB_;
  /** ratio per particle per sample
   *  size: particles * samples
   *  vector is over species
   *  each matrix row: particle col: sample
   */
  std::vector<Matrix<Value>> Psi_NM_;
  std::vector<Matrix<Value>> Phi_Psi_NB_, N_BB_;
  /** basis_values_ at each r of rsamples_ row: sample col: basis_value
   *  size: samples * basis_size
   */
  Matrix<Value> Phi_MB_;
  /** @} */

  /** @ingroup DensityIntegration only used for density integration
   *  @{
   */
  /// number of moves
  int nmoves_ = 0;
  /// number of accepted samples
  int naccepted_ = 0;
  /// running acceptance ratio over all samples
  Real acceptance_ratio_ = 0.0;
  int ind_dims_[OHMMS_DIM];
  bool warmed_up_ = false;
  /// }@

  Real metric_ = 1.0;

  // \todo is this state necessay, would it be better passed down the call stack?
  /// current position -- As long Positions are TinyVectors they are intialized to zero vectors
  Position rpcur_;
  /// current drift
  Position dpcur_;
  /// current density
  Real rhocur_ = 0.0;

public:
  /** Standard Constructor
   *  Call this to make a new OBDM this is what you should be calling
   */
  OneBodyDensityMatrices(OneBodyDensityMatricesInput&& obdmi,
                         const Lattice& lattice,
                         const SpeciesSet& species,
                         const WaveFunctionFactory& wf_factory,
                         ParticleSet& pset_target);

  /** Constructor used when spawing crowd clones
   *  needs to be public so std::make_unique can call it.
   *  Do not use directly unless you've really thought it through.
   */
  OneBodyDensityMatrices(const OneBodyDensityMatrices& obdm, DataLocality dl);

  std::unique_ptr<OperatorEstBase> spawnCrowdClone() const override;

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator& rng) override;

  void startBlock(int steps) override;

  /** create and tie OperatorEstimator's observable_helper hdf5 wrapper to stat.h5 file
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  void registerOperatorEstimator(hid_t gid) override;

private:
  /** Default copy constructor.
   *  Instances of this estimator is assume to be thread scope, i.e. never
   *  called by more than one thread at a time. note the OperatorEstBase copy constructor does
   *  not copy or even allocate data_
   */
  OneBodyDensityMatrices(const OneBodyDensityMatrices& obdm) = default;

  /** Unfortunate design RandomGenerator type aliasing and
   *  virtual inheritance requires this for testing.
   */
  template<class RNG_GEN>
  void implAccumulate(const RefVector<MCPWalker>& walkers,
                      const RefVector<ParticleSet>& psets,
                      const RefVector<TrialWaveFunction>& wfns,
                      RNG_GEN& rng);

  size_t calcFullDataSize(size_t basis_size, int num_species);
  //local functions
  void normalizeBasis(ParticleSet& pset_target);
  //  printing
  void report(const std::string& pad = "");
  template<class RNG_GEN>
  void evaluateMatrix(ParticleSet& pset_target, TrialWaveFunction& psi_target, const MCPWalker& walker, RNG_GEN& rng);
  //  sample generation
  /** Dispatch method to difference methods of generating samples.
   *  dispatch determined by Integrator.
   *  \param[in] weight       of this walker's samples
   *  \param[in] pset_target  will be returned to its initial state but is mutated.
   *  \param[in] rng          random generator. templated for testing without dependency on app level rng.
   *  \param[in] steps        If integrator_ = Integrator::DENSITY steps is a key parameter otherwise ignored.
   *                          when set to 0 it is reset to samples_ internally
   *  
   *  sideeffects:
   *   * samples_weights_ are set.
   *   * rsamples_ are set.
   *     for Density
   *      * update basis_values_, basis_gradients_, basis_laplacians_
   */
  // These functions deserve unit tests and likely should be pure functions.
  template<class RNG_GEN>
  void generateSamples(const Real weight, ParticleSet& pset_target, RNG_GEN& rng, int steps = 0);
  template<class RNG_GEN>
  void generateUniformGrid(RNG_GEN& rng);
  template<class RNG_GEN>
  void generateUniformSamples(RNG_GEN& rng);
  /** generate samples for density integration
   *  \param[in]   save          if false throw out the samples
   *  \param[in]   steps         actually the number of samples which are basically steps.
   *  \param[in]   rng           random generator
   *  \param[in]   pset_target   will be returned to its initial state but is mutated.
   *
   *  sideeffects:
   *   *
   */
  template<class RNG_GEN>
  void generateDensitySamples(bool save, int steps, RNG_GEN& rng, ParticleSet& pset_target);
  void generateSampleRatios(ParticleSet& pset_target,
                            TrialWaveFunction& psi_target,
                            std::vector<Matrix<Value>>& Psi_nm);
  /// produce a position difference vector from timestep
  template<class RNG_GEN>
  Position diffuse(const Real sqt, RNG_GEN& rng);
  /** calculate density based on r
   *  \param[in]      r       position
   *  \param[out]   dens      density
   *
   *  called by generateDensitySamples to get trial dens.
   *  also called by test_derivatives.
   *  sideeffects:
   *    * updateBasis is called
   */
  void calcDensity(const Position& r, Real& dens, ParticleSet& pset_target);
  /** calculate density and drift bashed on r
   *  \param[in]      r       position
   *  \param[out]   dens      density
   *  \param[out]   drift     drift
   *
   *  called by warmupSamples to get an initial drift and density.
   *  called by generateDensitySamples to get trial drift and trial dens.
   *  also called by test_derivatives.
   *  sideeffects:
   *    * updateBasisD012 is called
   */
  void calcDensityDrift(const Position& r, Real& dens, Position& drift, ParticleSet& pset_target);
  //  basis & wavefunction ratio matrix construction

  /** set Phi_mp to basis vaules per sample
   *  sideeffects:
   *    * updates basis_values_ to last rsample
   */
  void generateSampleBasis(Matrix<Value>& Phi_mb, ParticleSet& pset_target, TrialWaveFunction& psi_target);
  /** set phi_nb to basis values per target particleset particle
   *  sideeffects:
   *    * updates basis_values_ to last rsample
   */
  void generateParticleBasis(ParticleSet& pset_target, std::vector<Matrix<Value>>& phi_nb);

  //  basis set updates
  void updateBasis(const Position& r, ParticleSet& pset_target);
  /** evaluates vgl on basis_functions_ for r
   *  sideeffects:
   *    * sets basis_values_, basis_gradients_, basis_laplacians_
   *      all are normalized by basis norms_
   */
  void updateBasisD012(const Position& r, ParticleSet& pset_target);
  /** does some warmup sampling i.e. samples but throws away the results
   *  Only when integrator_ = Integrator::DENSITY
   *  sets rpcur_ initial rpcur + one diffusion step
   *  sets initial rhocur_ and dpcur_
   *  Then calls generateSamples with number of input warmup samples.
   */
  template<typename RAN_GEN>
  void warmupSampling(ParticleSet& pset_target, RAN_GEN& rng);

  struct OneBodyDensityMatrixTimers
  {
    NewTimer& eval_timer;
    NewTimer& gen_samples_timer;
    NewTimer& gen_sample_basis_timer;
    NewTimer& gen_sample_ratios_timer;
    NewTimer& gen_particle_basis_timer;
    NewTimer& matrix_products_timer;
    NewTimer& accumulate_timer;
    OneBodyDensityMatrixTimers(const std::string& prefix)
        : eval_timer(*timer_manager.createTimer(prefix + "Eval", timer_level_fine)),
          gen_samples_timer(*timer_manager.createTimer(prefix + "GenSamples", timer_level_fine)),
          gen_sample_basis_timer(*timer_manager.createTimer(prefix + "GenSampleBasis", timer_level_fine)),
          gen_sample_ratios_timer(*timer_manager.createTimer(prefix + "GenSampleRatios", timer_level_fine)),
          gen_particle_basis_timer(*timer_manager.createTimer(prefix + "GenParticleBasis", timer_level_fine)),
          matrix_products_timer(*timer_manager.createTimer(prefix + "MatrixProducts", timer_level_fine)),
          accumulate_timer(*timer_manager.createTimer(prefix + "Accumulate", timer_level_fine))
    {}
  };

  OneBodyDensityMatrixTimers timers_;

public:
  template<typename T>
  friend class testing::OneBodyDensityMatricesTests;
};

extern template void OneBodyDensityMatrices::generateSamples<RandomGenerator>(Real weight,
                                                                              ParticleSet& pset_target,
                                                                              RandomGenerator& rng,
                                                                              int steps);
extern template void OneBodyDensityMatrices::evaluateMatrix<RandomGenerator>(ParticleSet& pset_target,
                                                                             TrialWaveFunction& psi_target,
                                                                             const MCPWalker& walker,
                                                                             RandomGenerator& rng);
extern template void OneBodyDensityMatrices::implAccumulate<RandomGenerator>(const RefVector<MCPWalker>& walkers,
                                                                             const RefVector<ParticleSet>& psets,
                                                                             const RefVector<TrialWaveFunction>& wfns,
                                                                             RandomGenerator& rng);
#if defined(USE_FAKE_RNG) || defined(QMC_RNG_BOOST)
extern template void OneBodyDensityMatrices::generateSamples<StdRandom<double>>(Real weight,
                                                                                ParticleSet& pset_target,
                                                                                StdRandom<double>& rng,
                                                                                int steps);
extern template void OneBodyDensityMatrices::evaluateMatrix<StdRandom<double>>(ParticleSet& pset_target,
                                                                               TrialWaveFunction& psi_target,
                                                                               const MCPWalker& walker,
                                                                               StdRandom<double>& rng);
extern template void OneBodyDensityMatrices::implAccumulate<StdRandom<double>>(const RefVector<MCPWalker>& walkers,
                                                                               const RefVector<ParticleSet>& psets,
                                                                               const RefVector<TrialWaveFunction>& wfns,
                                                                               StdRandom<double>& rng);
#endif

} // namespace qmcplusplus

#endif
