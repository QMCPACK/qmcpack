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
  /// actual WaveFunctionFactory instance must be owned by the same or enclosing scope
  const WaveFunctionFactory& wf_factory_;

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

  //data members \todo remove unecessary state
  CompositeSPOSet basis_functions;
  Vector<Value> basis_values;
  Vector<Value> basis_norms;
  Vector<Grad> basis_gradients;
  Vector<Value> basis_laplacians;
  std::vector<Position> rsamples;
  Vector<Real> sample_weights;
  std::vector<Value> psi_ratios;
  Real dens;
  Position drift;
  int nindex;
  bool volume_normed;
  int basis_size;
  std::vector<int> species_size;
  std::vector<std::string> species_name;
  std::vector<Matrix<Value>> Phi_NB, Psi_NM, Phi_Psi_NB, N_BB, E_BB;
  Matrix<Value> Phi_MB;
  bool check_overlap;
  bool check_derivatives;

// \todo this seems like developement clutter remove?
//#define DMCHECK
#ifdef DMCHECK
  std::vector<Vector<Value>*> E_Ntmp;
  std::vector<Matrix<Value>*> Phi_NBtmp, Psi_NMtmp, Phi_Psi_NBtmp, N_BBtmp, E_BBtmp;
  Matrix_t Phi_MBtmp;
#endif

  /** @ingroup DensityIntegration only used for density integration
   *  @{
   */
  /// number of moves
  int nmoves = 0;
  /// number of accepted samples
  int naccepted = 0;
  /// running acceptance ratio over all samples
  Real acceptance_ratio = 0.0;
  bool write_acceptance_ratio;
  int ind_dims_[OHMMS_DIM];
  /// }@

  Real metric_;

  /// current position
  Position rpcur;
  /// current drift
  Position dpcur;
  /// current density
  Real rhocur;

public:
  /** Standard Constructor
   *  If you are making a new OBDM this is what you should be calling
   */
  OneBodyDensityMatrices(OneBodyDensityMatricesInput&& obdmi,
                         const Lattice& lattice,
                         const SpeciesSet& species,
                         const WaveFunctionFactory& wf_factory);

  /** copy constructor delegates to standard constructor
   *  This results in a copy construct and move of OneBodyDensityMatricesInput
   *  But for the OBDM itself its as if it went through the standard construction.
   *  This could be optimized.
   */
  OneBodyDensityMatrices(const OneBodyDensityMatrices& obdm);
  ~OneBodyDensityMatrices() override;

  OneBodyDensityMatrices* clone() override;

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
  //local functions
  void normalize(ParticleSet& pset_target);
  //  printing
  void report(const std::string& pad = "");
  //  sample generation
  template<class RNG_GEN>
  void generateSamples(Real weight, ParticleSet& pset_target, RNG_GEN& rng, int steps = 0);
  template<class RNG_GEN>
  void generate_uniform_grid(RNG_GEN& rng);
  template<class RNG_GEN>
  void generate_uniform_samples(RNG_GEN& rng);
  template<class RNG_GEN>
  void generate_density_samples(bool save, int steps, RNG_GEN& rng, ParticleSet& pset_target);
  /// produce a position difference vector from timestep
  template<class RNG_GEN>
  Position diffusion(const Real sqt, RNG_GEN& rng);
  void density_only(const Position& r, Real& dens, ParticleSet& pset_target);
  void density_drift(const Position& r, Real& dens, Position& drift, ParticleSet& pset_target);
  //  basis & wavefunction ratio matrix construction
  void get_energies(std::vector<Vector<Value>*>& E_n);
  void generate_sample_basis(Matrix<Value>& Phi_mb, ParticleSet& pset_target, TrialWaveFunction& psi_target);
  void generate_sample_ratios(std::vector<Matrix<Value>>& Psi_nm,
                              ParticleSet& pset_target,
                              TrialWaveFunction& psi_target);
  void generate_particle_basis(ParticleSet& P, std::vector<Matrix<Value>>& Phi_nb, ParticleSet& pset_target);
  //  basis set updates
  void update_basis(const Position& r, ParticleSet& pset_target);
  void update_basis_d012(const Position& r, ParticleSet& pset_target);

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

/** OneBodyDensityMatrices factory function
 */
UPtr<OneBodyDensityMatrices> createOneBodyDensityMatrices(OneBodyDensityMatricesInput&& obdmi,
                                                          const PtclOnLatticeTraits::ParticleLayout_t lattice,
                                                          const SpeciesSet& species);


extern template void OneBodyDensityMatrices::generateSamples<RandomGenerator_t>(Real weight,
                                                                                ParticleSet& pset_target,
                                                                                RandomGenerator_t& rng,
                                                                                int steps);
extern template void OneBodyDensityMatrices::generateSamples<StdRandom<double>>(Real weight,
                                                                                ParticleSet& pset_target,
                                                                                StdRandom<double>& rng,
                                                                                int steps);
} // namespace qmcplusplus

#endif
