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

namespace qmcplusplus
{

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
  Sampling sampling_;
  //data members
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
  int samples_;
  int nparticles;
  int nspecies;
  std::vector<int> species_size;
  std::vector<std::string> species_name;
  std::vector<Matrix<Value>*> Phi_NB, Psi_NM, Phi_Psi_NB, N_BB, E_BB;
  Matrix<Value> Phi_MB;
  bool check_overlap;
  bool check_derivatives;

//#define DMCHECK
#ifdef DMCHECK
  std::vector<Vector<Value>*> E_Ntmp;
  std::vector<Matrix<Value>*> Phi_NBtmp, Psi_NMtmp, Phi_Psi_NBtmp, N_BBtmp, E_BBtmp;
  Matrix_t Phi_MBtmp;
#endif

  Position center_;
  Position rcorner_;
  Real volume_;
  bool periodic_;
  int nmoves;
  int naccepted;
  Real acceptance_ratio;
  bool write_acceptance_ratio;
  bool write_rstats;

  int ind_dims_[OHMMS_DIM];
  Real metric;

  Position rpcur;
  Position dpcur;
  Real rhocur;

  RandomGenerator_t* uniform_random;

public:
  /** Standard Constructor
   *  If you are making a new OBDM this is what you should be calling
   */
  OneBodyDensityMatrices(OneBodyDensityMatricesInput&& obdmi, const Lattice& lattice, const SpeciesSet& species);

  /** actual copy constructor
   */
  OneBodyDensityMatrices(const OneBodyDensityMatrices& obdm);
  ~OneBodyDensityMatrices() override;

  OneBodyDensityMatrices* clone() override;

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator_t& rng) override;

  void normalize(Real invToWgt) override;

  void setRandomGenerator(RandomGenerator_t* rng);

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
  void normalize();
  //  printing
  void report(const std::string& pad = "");
  //  sample generation
  void generate_samples(Real weight, ParticleSet& pset_target, int steps = 0);
  void generate_uniform_grid(RandomGenerator_t& rng);
  void generate_uniform_samples(RandomGenerator_t& rng);
  void generate_density_samples(bool save, int steps, RandomGenerator_t& rng, ParticleSet& pset_target);
  void diffusion(Real sqt, Position& diff);
  void density_only(const Position& r, Real& dens, ParticleSet& pset_target);
  void density_drift(const Position& r, Real& dens, Position& drift, ParticleSet& pset_target);
  //  basis & wavefunction ratio matrix construction
  void get_energies(std::vector<Vector<Value>*>& E_n);
  void generate_sample_basis(Matrix<Value>& Phi_mb, ParticleSet& pset_target, TrialWaveFunction& psi_target);
  void generate_sample_ratios(std::vector<Matrix<Value>*> Psi_nm,
                              ParticleSet& pset_target,
                              TrialWaveFunction& psi_target);
  void generate_particle_basis(ParticleSet& P, std::vector<Matrix<Value>*>& Phi_nb, ParticleSet& pset_target);
  //  basis set updates
  void update_basis(const Position& r, ParticleSet& pset_target);
  void update_basis_d012(const Position& r, ParticleSet& pset_target);
  void normalize(ParticleSet& pset_target);

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
};

/** OneBodyDensityMatrices factory function
 */
UPtr<OneBodyDensityMatrices> createOneBodyDensityMatrices(OneBodyDensityMatricesInput&& obdmi,
                                                          const PtclOnLatticeTraits::ParticleLayout_t lattice,
                                                          const SpeciesSet& species);


} // namespace qmcplusplus

#endif
