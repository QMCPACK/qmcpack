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

//template<typename VALUE = QMCTraits::ValueType, typename VALUE_FP = QMCTraits::FullPrecValueType>
class OneBodyDensityMatrices : public OperatorEstBase
{
public:
  using Value         = QMCTraits::ValueType;
  using FullPrecValue = QMCTraits::FullPrecValueType;
  using Real          = RealAlias<Value>;
  using FullPrecReal  = RealAlias<FullPrecValue>;
  using Grad          = TinyVector<Value, OHMMS_DIM>;
  using Lattice       = PtclOnLatticeTraits::ParticleLayout_t;
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
  /** WaveFunctionFactory reference to allow delegation of the copy constructor
   *  \todo remove after copy constructor that directly shares or copys basis_set_ is done
   */
  const WaveFunctionFactory& wf_factory_;
  /** target particleset  reference to allow delegation of the copy constructor
   *  \todo remove after copy constructor that directly shares or copys basis_set_ is done
   */
  ParticleSet& very_temp_pset_;

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
  std::vector<Value> psi_ratios_;
  int basis_size_;
  std::vector<int> species_sizes_;
  std::vector<std::string> species_names_;
  std::vector<Matrix<Value>> Phi_NB_, Psi_NM_, Phi_Psi_NB_, N_BB_;
  Matrix<Value> Phi_MB_;

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

  Real metric_;

  // \todo is this state necessay, would it be better passed down the call stack?
  /// current position
  Position rpcur_;
  /// current drift
  Position dpcur_;
  /// current density
  Real rhocur_;

public:
  /** Standard Constructor
   *  If you are making a new OBDM this is what you should be calling
   */
  OneBodyDensityMatrices(OneBodyDensityMatricesInput&& obdmi,
                         const Lattice& lattice,
                         const SpeciesSet& species,
                         const WaveFunctionFactory& wf_factory,
                         ParticleSet& pset_target,
                         const DataLocality dl = DataLocality::crowd);

  /** copy constructor delegates to standard constructor
   *  This results in a copy construct and move of OneBodyDensityMatricesInput
   *  But for the OBDM itself its as if it went through the standard construction.
   *  This will be replaced within a few PR's by an optimized copy constructor.
   */
  OneBodyDensityMatrices(const OneBodyDensityMatrices& obdm);
  ~OneBodyDensityMatrices() override;

  std::unique_ptr<OperatorEstBase> clone() const override;

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator_t& rng) override;

  void normalize(Real invToWgt) override;
  void startBlock(int steps) override;

  /** create and tie OperatorEstimator's observable_helper hdf5 wrapper to stat.h5 file
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  void registerOperatorEstimator(hid_t gid) override {}

private:
  size_t calcFullDataSize(size_t basis_size, int num_species);
  //local functions
  void normalize(ParticleSet& pset_target);
  //  printing
  void report(const std::string& pad = "");
  template<class RNG_GEN>
  void evaluateMatrix(ParticleSet& pset_target, TrialWaveFunction& psi_target, const MCPWalker& walker, RNG_GEN& rng);
  //  sample generation
  template<class RNG_GEN>
  void generateSamples(Real weight, ParticleSet& pset_target, RNG_GEN& rng, int steps = 0);
  // These functions deserve unit tests and likely should be pure functions.
  template<class RNG_GEN>
  void generateUniformGrid(RNG_GEN& rng);
  template<class RNG_GEN>
  void generateUniformSamples(RNG_GEN& rng);
  template<class RNG_GEN>
  void generateDensitySamples(bool save, int steps, RNG_GEN& rng, ParticleSet& pset_target);
  void generateSampleRatios(ParticleSet& pset_target,
                            TrialWaveFunction& psi_target,
                            std::vector<Matrix<Value>>& Psi_nm);
  /// produce a position difference vector from timestep
  template<class RNG_GEN>
  Position diffuse(const Real sqt, RNG_GEN& rng);
  void calcDensity(const Position& r, Real& dens, ParticleSet& pset_target);
  void calcDensityDrift(const Position& r, Real& dens, Position& drift, ParticleSet& pset_target);
  //  basis & wavefunction ratio matrix construction
  void generateSampleBasis(Matrix<Value>& Phi_mb, ParticleSet& pset_target, TrialWaveFunction& psi_target);
  void generateParticleBasis(ParticleSet& pset_target, std::vector<Matrix<Value>>& phi_nb);

  //  basis set updates
  void updateBasis(const Position& r, ParticleSet& pset_target);
  void updateBasisD012(const Position& r, ParticleSet& pset_target);
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

extern template void OneBodyDensityMatrices::generateSamples<RandomGenerator_t>(Real weight,
                                                                                ParticleSet& pset_target,
                                                                                RandomGenerator_t& rng,
                                                                                int steps);
extern template void OneBodyDensityMatrices::generateSamples<StdRandom<double>>(Real weight,
                                                                                ParticleSet& pset_target,
                                                                                StdRandom<double>& rng,
                                                                                int steps);
extern template void OneBodyDensityMatrices::evaluateMatrix<RandomGenerator_t>(ParticleSet& pset_target,
                                                                               TrialWaveFunction& psi_target,
                                                                               const MCPWalker& walker,
                                                                               RandomGenerator_t& rng);
extern template void OneBodyDensityMatrices::evaluateMatrix<StdRandom<double>>(ParticleSet& pset_target,
                                                                               TrialWaveFunction& psi_target,
                                                                               const MCPWalker& walker,
                                                                               StdRandom<double>& rng);

} // namespace qmcplusplus

#endif
