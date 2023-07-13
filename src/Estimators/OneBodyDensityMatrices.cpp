//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCHamiltonians/DensityMatrices1B.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "OneBodyDensityMatrices.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Numerics/MatrixOperators.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/string_utils.h"
#include "type_traits/complex_help.hpp"
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
using MatrixOperators::diag_product;
using MatrixOperators::product;
using MatrixOperators::product_AtB;

OneBodyDensityMatrices::OneBodyDensityMatrices(OneBodyDensityMatricesInput&& obdmi,
                                               const Lattice& lattice,
                                               const SpeciesSet& species,
                                               const SPOMap& spomap,
                                               const ParticleSet& pset_target)
    : OperatorEstBase(DataLocality::crowd),
      input_(obdmi),
      lattice_(lattice),
      species_(species),
      basis_functions_("OneBodyDensityMatrices::basis"),
      timers_("OneBodyDensityMatrix")
{
  my_name_ = "OneBodyDensityMatrices";
  lattice_.reset();

  if (input_.get_corner_defined())
  {
    rcorner_ = input_.get_corner();
    center_  = rcorner_ + input_.get_scale() * lattice_.Center;
  }
  else
  {
    if (input_.get_center_defined())
      center_ = input_.get_center();
    else
      center_ = lattice_.Center;
    rcorner_ = center_ - input_.get_scale() * lattice_.Center;
  }

  volume_   = lattice_.Volume * std::exp(OHMMS_DIM * std::log(input_.get_scale()));
  periodic_ = lattice_.SuperCellEnum != SUPERCELL_OPEN;

  // Here we discover sampling is derived (this may belong in input class)
  switch (input_.get_integrator())
  {
  case Integrator::UNIFORM_GRID:
    sampling_    = Sampling::VOLUME_BASED;
    samples_     = pow(input_.get_points(), OHMMS_DIM);
    metric_      = volume_ / samples_;
    ind_dims_[0] = pow(input_.get_points(), OHMMS_DIM - 1);
    for (int d = 1; d < OHMMS_DIM; ++d)
      ind_dims_[d] = ind_dims_[d - 1] / input_.get_points();
    break;
  case Integrator::UNIFORM:
    sampling_ = Sampling::VOLUME_BASED;
    samples_  = input_.get_samples();
    metric_   = volume_ / samples_;
    break;
  case Integrator::DENSITY:
    sampling_ = Sampling::METROPOLIS;
    samples_  = input_.get_samples();
    metric_   = 1.0 / samples_;
    break;
  }
  rsamples_.resize(samples_);

  // get the sposets that form the basis
  auto& sposets = input_.get_basis_sets();

  for (int i = 0; i < sposets.size(); ++i)
  {
    auto spo_it = spomap.find(sposets[i]);
    if (spo_it == spomap.end())
      throw UniformCommunicateError("OneBodyDensityMatrices::OneBodyDensityMatrices sposet " + sposets[i] +
                                    " does not exist.");
    basis_functions_.add(spo_it->second->makeClone());
  }
  basis_size_ = basis_functions_.size();

  if (basis_size_ < 1)
    throw UniformCommunicateError("OneBodyDensityMatrices::OneBodyDensityMatrices basis_size must be greater than one");

  int nspecies = species.size();
  if (!species_.hasAttribute("membersize"))
    throw UniformCommunicateError("OneBodyDensityMatrices::OneBodyDensityMatrices error: Species set does not have the "
                                  "required attribute 'membersize'");
  int isize = species.getAttribute("membersize");
  // We have the count per species at least a fundamental as the  total particles.
  // the sum of species membersize and total particles should be an invariant.
  int nparticles = 0;
  for (int s = 0; s < nspecies; ++s)
    nparticles += species(isize, s);
  for (int s = 0; s < nspecies; ++s)
    species_sizes_.push_back(species(isize, s));
  for (int s = 0; s < nspecies; ++s)
    species_names_.push_back(species.speciesName[s]);

  basis_values_.resize(basis_size_);
  basis_norms_.resize(basis_size_);

  Real bn_standard = 1.0;
  if (input_.get_volume_normalized())
    bn_standard = 1.0 / std::sqrt(volume_);
  for (int i = 0; i < basis_size_; ++i)
    basis_norms_[i] = bn_standard;

  rsamples_.resize(samples_);
  samples_weights_.resize(samples_);
  psi_ratios_.resize(nparticles);

  if (input_.get_evaluator() == Evaluator::MATRIX)
  {
    Phi_MB_.resize(samples_, basis_size_);
    Phi_NB_.reserve(nspecies);
    Psi_NM_.reserve(nspecies);
    Phi_Psi_NB_.reserve(nspecies);
    N_BB_.reserve(nspecies);
    for (int s = 0; s < nspecies; ++s)
    {
      int specs_size = species_sizes_[s];
      Phi_NB_.emplace_back(specs_size, basis_size_);
      Psi_NM_.emplace_back(specs_size, samples_);
      Phi_Psi_NB_.emplace_back(specs_size, basis_size_);
      N_BB_.emplace_back(basis_size_, basis_size_);
    }
  }

  if (sampling_ == Sampling::METROPOLIS)
  {
    basis_gradients_.resize(basis_size_);
    basis_laplacians_.resize(basis_size_);
  }

  // so if the input is not normalized, normalize it.
  // with respect to what?
  if (!input_.get_normalized())
  {
    //Since the following is not a const method we copy particle set
    ParticleSet pset_temp(pset_target);
    normalizeBasis(pset_temp);
  }

  data_.resize(calcFullDataSize(basis_size_, species_.size()), 0.0);
}

OneBodyDensityMatrices::OneBodyDensityMatrices(const OneBodyDensityMatrices& obdm, DataLocality dl)
    : OneBodyDensityMatrices(obdm)
{
  data_locality_ = dl;
}

std::unique_ptr<OperatorEstBase> OneBodyDensityMatrices::spawnCrowdClone() const
{
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;

  if (data_locality_ == DataLocality::rank)
  {
    // This is just a stub until a memory saving optimization is deemed necessary
    spawn_data_locality = DataLocality::queue;
    data_size           = 0;
    throw std::runtime_error("There is no memory savings implementation for OneBodyDensityMatrices");
  }

  auto spawn = std::make_unique<OneBodyDensityMatrices>(*this, spawn_data_locality);
  spawn->get_data().resize(data_size, 0.0);
  return spawn;
}

size_t OneBodyDensityMatrices::calcFullDataSize(const size_t basis_size, const int nspecies)
{
  if constexpr (IsComplex_t<Value>::value)
    return 2 * basis_size * basis_size * nspecies;
  else
    return basis_size * basis_size * nspecies;
}

void OneBodyDensityMatrices::startBlock(int steps) {}

void OneBodyDensityMatrices::generateSamples(const Real weight,
                                             ParticleSet& pset_target,
                                             RandomBase<FullPrecReal>& rng,
                                             int steps)
{
  ScopedTimer local_timer(timers_.gen_samples_timer);

  // Steps will always be 0 unless these are samples for warmup which is only for metropolis
  // This is not a clear way to write this
  // \todo rewrite to make algorithm more clears
  bool save = false;
  if (steps == 0)
  {
    save  = true;
    steps = samples_;
  }

  switch (input_.get_integrator())
  {
  case Integrator::UNIFORM_GRID:
    generateUniformGrid(rng);
    break;
  case Integrator::UNIFORM:
    generateUniformSamples(rng);
    break;
  case Integrator::DENSITY: {
    generateDensitySamples(save, steps, rng, pset_target);
  }
  }

  if (save)
  {
    if (sampling_ == Sampling::METROPOLIS)
      samples_weights_ *= weight;
    else
    {
      std::fill(samples_weights_.begin(), samples_weights_.end(), weight);
    }
  }

  // optional check
  if (input_.get_rstats() && omp_get_thread_num() == 0)
  {
    Position rmin  = std::numeric_limits<Real>::max();
    Position rmax  = -std::numeric_limits<Real>::max();
    Position rmean = 0.0;
    Position rstd  = 0.0;
    for (int s = 0; s < rsamples_.size(); ++s)
      for (int d = 0; d < OHMMS_DIM; ++d)
      {
        Real rd = rsamples_[s][d];
        rmin[d] = std::min(rmin[d], rd);
        rmax[d] = std::max(rmax[d], rd);
        rmean[d] += rd;
        rstd[d] += rd * rd;
      }
    rmean /= rsamples_.size();
    rstd /= rsamples_.size();
    for (int d = 0; d < OHMMS_DIM; ++d)
      rstd[d] = std::sqrt(rstd[d] - rmean[d] * rmean[d]);
    app_log() << "\nrsamples properties:" << std::endl;
    app_log() << "  rmin  = " << rmin << std::endl;
    app_log() << "  rmax  = " << rmax << std::endl;
    app_log() << "  rmean = " << rmean << std::endl;
    app_log() << "  rstd  = " << rstd << std::endl;
  }
}

inline void OneBodyDensityMatrices::generateUniformGrid(RandomBase<FullPrecReal>& rng)
{
  Position rp;
  Position ushift = 0.0;
  Real du         = input_.get_scale() / input_.get_points();
  for (int d = 0; d < OHMMS_DIM; ++d)
    ushift[d] += rng() * du;
  for (int s = 0; s < samples_; ++s)
  {
    int nrem = s;
    for (int d = 0; d < OHMMS_DIM - 1; ++d)
    {
      int ind = nrem / ind_dims_[d];
      rp[d]   = ind * du + ushift[d];
      nrem -= ind * ind_dims_[d];
    }
    rp[OHMMS_DIM - 1] = nrem * du + ushift[OHMMS_DIM - 1];
    rsamples_[s]      = lattice_.toCart(rp) + rcorner_;
  }
}

inline void OneBodyDensityMatrices::generateUniformSamples(RandomBase<FullPrecReal>& rng)
{
  Position rp;
  for (int s = 0; s < samples_; ++s)
  {
    for (int d = 0; d < OHMMS_DIM; ++d)
      rp[d] = input_.get_scale() * rng();
    rsamples_[s] = lattice_.toCart(rp) + rcorner_;
  }
}

inline void OneBodyDensityMatrices::generateDensitySamples(bool save,
                                                           int steps,
                                                           RandomBase<FullPrecReal>& rng,
                                                           ParticleSet& pset_target)
{
  const auto timestep = input_.get_timestep();
  Real sqt            = std::sqrt(timestep);
  Real ot             = 1.0 / timestep;
  Position r          = rpcur_;  //current position
  Position d          = dpcur_;  //current drift
  Real rho            = rhocur_; //current density
  for (int s = 0; s < steps; ++s)
  {
    nmoves_++;
    Position rp;                       // trial pos
    Position dp;                       // trial drift
    Position ds;                       // drift sum
    Real rhop;                         // trial density
    Real ratio;                        // dens ratio
    Real Pacc;                         // acc prob
    Position diff = diffuse(sqt, rng); // get diffusion
    if (input_.get_use_drift())
    {
      rp = r + diff + d;                                                  //update trial position
      calcDensityDrift(rp, rhop, dp, pset_target);                        //get trial drift and density
      ratio = rhop / rho;                                                 //density ratio
      ds    = dp + d;                                                     //drift sum
      Pacc  = ratio * std::exp(-ot * (dot(diff, ds) + .5 * dot(ds, ds))); //acceptance probability
    }
    else
    {
      rp = r + diff;                      //update trial position
      calcDensity(rp, rhop, pset_target); //get trial density
      ratio = rhop / rho;                 //density ratio
      Pacc  = ratio;                      //acceptance probability
    }
    if (rng() < Pacc)
    { //accept move
      r   = rp;
      d   = dp;
      rho = rhop;
      naccepted_++;
    }
    if (save)
    {
      rsamples_[s]        = r;
      samples_weights_[s] = 1.0 / rho;
    }
  }
  acceptance_ratio_ = Real(naccepted_) / nmoves_;

  if (input_.get_write_acceptance_ratio() && omp_get_thread_num() == 0)
    app_log() << "dm1b  acceptance_ratio = " << acceptance_ratio_ << std::endl;

  rpcur_  = r;
  dpcur_  = d;
  rhocur_ = rho;
}

OneBodyDensityMatrices::Position OneBodyDensityMatrices::diffuse(const Real sqt, RandomBase<FullPrecReal>& rng)
{
  Position diff;
  assignGaussRand(&diff[0], OHMMS_DIM, rng);
  diff *= sqt;
  return diff;
}


inline void OneBodyDensityMatrices::calcDensity(const Position& r, Real& dens, ParticleSet& pset_target)
{
  updateBasis(r, pset_target);
  dens = 0.0;
  for (int i = 0; i < basis_size_; ++i)
  {
    Value b = basis_values_[i];
    dens += std::abs(qmcplusplus::conj(b) * b);
  }
  dens /= basis_size_;
}


void OneBodyDensityMatrices::calcDensityDrift(const Position& r, Real& dens, Position& drift, ParticleSet& pset_target)
{
  updateBasisD012(r, pset_target);
  dens  = 0.0;
  drift = 0.0;
  for (int i = 0; i < basis_size_; ++i)
  {
    const Grad& bg = basis_gradients_[i];
    Value b        = basis_values_[i];
    Value bc       = qmcplusplus::conj(b);
    dens += std::abs(bc * b);
    for (int d = 0; d < OHMMS_DIM; ++d)
      drift[d] += std::real(bc * bg[d]);
  }
  drift *= input_.get_timestep() / dens;
  dens /= basis_size_;
}

void OneBodyDensityMatrices::accumulate(const RefVector<MCPWalker>& walkers,
                                        const RefVector<ParticleSet>& psets,
                                        const RefVector<TrialWaveFunction>& wfns,
                                        RandomBase<FullPrecReal>& rng)
{
  implAccumulate(walkers, psets, wfns, rng);
}

void OneBodyDensityMatrices::implAccumulate(const RefVector<MCPWalker>& walkers,
                                            const RefVector<ParticleSet>& psets,
                                            const RefVector<TrialWaveFunction>& wfns,
                                            RandomBase<FullPrecReal>& rng)
{
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    walkers_weight_ += walkers[iw].get().Weight;
    evaluateMatrix(psets[iw], wfns[iw], walkers[iw], rng);
  }
}

void OneBodyDensityMatrices::evaluateMatrix(ParticleSet& pset_target,
                                            TrialWaveFunction& psi_target,
                                            const MCPWalker& walker,
                                            RandomBase<FullPrecReal>& rng)
{
  //perform warmup sampling the first time
  warmupSampling(pset_target, rng);
  // get weight and single particle energy trace data
  Real weight = walker.Weight * metric_;

  // compute sample positions (monte carlo or deterministic)
  generateSamples(weight, pset_target, rng);
  // compute basis and wavefunction ratio values in matrix form
  generateSampleBasis(Phi_MB_, pset_target, psi_target);  // basis           : samples   x basis_size
  generateSampleRatios(pset_target, psi_target, Psi_NM_); // conj(Psi ratio) : particles x samples
  generateParticleBasis(pset_target, Phi_NB_);            // conj(basis)     : particles x basis_size

  // \todo separate testable and optimizable block, should be function
  // perform integration via matrix products
  {
    ScopedTimer local_timer(timers_.matrix_products_timer);
    for (int s = 0; s < species_.size(); ++s)
    {
      Matrix<Value>& Psi_nm     = Psi_NM_[s];
      Matrix<Value>& Phi_Psi_nb = Phi_Psi_NB_[s];
      Matrix<Value>& Phi_nb     = Phi_NB_[s];
      diag_product(Psi_nm, samples_weights_, Psi_nm);
      product(Psi_nm, Phi_MB_, Phi_Psi_nb);      // ratio*basis : particles x basis_size
      product_AtB(Phi_nb, Phi_Psi_nb, N_BB_[s]); // conj(basis)^T*ratio*basis : basis_size^2
    }
  }
  // accumulate data for this walker
  {
    ScopedTimer local_timer(timers_.accumulate_timer);
    const int basis_size_sq = basis_size_ * basis_size_;
    int ij                  = 0;
    for (int s = 0; s < species_.size(); ++s)
    {
      //int ij=nindex; // for testing
      const Matrix<Value>& NDM = N_BB_[s];
      for (int n = 0; n < basis_size_sq; ++n)
      {
        Value val = NDM(n);
        data_[ij] += real(val);
        ij++;
#if defined(QMC_COMPLEX)
        data_[ij] += imag(val);
        ij++;
#endif
      }
    }
  }
}

void OneBodyDensityMatrices::generateParticleBasis(ParticleSet& pset_target, std::vector<Matrix<Value>>& phi_nb)
{
  ScopedTimer local_timer(timers_.gen_particle_basis_timer);
  int p = 0;
  for (int s = 0; s < species_.size(); ++s)
  {
    int nb              = 0;
    Matrix<Value>& P_nb = phi_nb[s];
    for (int n = 0; n < species_sizes_[s]; ++n, ++p)
    {
      updateBasis(pset_target.R[p], pset_target);
      for (int b = 0; b < basis_size_; ++b, ++nb)
        P_nb(nb) = qmcplusplus::conj(basis_values_[b]);
    }
  }
}

void OneBodyDensityMatrices::generateSampleBasis(Matrix<Value>& Phi_mb,
                                                 ParticleSet& pset_target,
                                                 TrialWaveFunction& psi_target)
{
  ScopedTimer local_timer(timers_.gen_sample_basis_timer);
  int mb = 0;
  for (int m = 0; m < samples_; ++m)
  {
    updateBasis(rsamples_[m], pset_target);
    for (int b = 0; b < basis_size_; ++b, ++mb)
      Phi_mb(mb) = basis_values_[b];
  }
}

void OneBodyDensityMatrices::generateSampleRatios(ParticleSet& pset_target,
                                                  TrialWaveFunction& psi_target,
                                                  std::vector<Matrix<Value>>& psi_nm)
{
  ScopedTimer local_timer(timers_.gen_sample_ratios_timer);
  for (int m = 0; m < samples_; ++m)
  {
    // get N ratios for the current sample point
    pset_target.makeVirtualMoves(rsamples_[m]);
    psi_target.evaluateRatiosAlltoOne(pset_target, psi_ratios_);

    // collect ratios into per-species matrices
    int p = 0;
    for (int s = 0; s < species_.size(); ++s)
    {
      Matrix<Value>& P_nm = psi_nm[s];
      for (int n = 0; n < species_sizes_[s]; ++n, ++p)
      {
        P_nm(n, m) = qmcplusplus::conj(psi_ratios_[p]);
      }
    }
  }
}

inline void OneBodyDensityMatrices::updateBasis(const Position& r, ParticleSet& pset_target)
{
  // This is ridiculous in the case of splines, still necessary for hybrid/LCAO
  pset_target.makeMove(0, r - pset_target.R[0]);
  basis_functions_.evaluateValue(pset_target, 0, basis_values_);
  pset_target.rejectMove(0);
  for (int i = 0; i < basis_size_; ++i)
    basis_values_[i] *= basis_norms_[i];
}


inline void OneBodyDensityMatrices::updateBasisD012(const Position& r, ParticleSet& pset_target)
{
  pset_target.makeMove(0, r - pset_target.R[0]);
  basis_functions_.evaluateVGL(pset_target, 0, basis_values_, basis_gradients_, basis_laplacians_);
  pset_target.rejectMove(0);
  for (int i = 0; i < basis_size_; ++i)
    basis_values_[i] *= basis_norms_[i];
  for (int i = 0; i < basis_size_; ++i)
    basis_gradients_[i] *= basis_norms_[i];
  for (int i = 0; i < basis_size_; ++i)
    basis_laplacians_[i] *= basis_norms_[i];
}

void OneBodyDensityMatrices::warmupSampling(ParticleSet& pset_target, RandomBase<FullPrecReal>& rng)
{
  if (sampling_ == Sampling::METROPOLIS)
  {
    if (!warmed_up_)
    {
      rpcur_ = diffuse(std::sqrt(input_.get_timestep()), rng);
      rpcur_ += center_;
      calcDensityDrift(rpcur_, rhocur_, dpcur_, pset_target);
    }
    generateSamples(1.0, pset_target, rng, input_.get_warmup_samples());
    warmed_up_ = true;
  }
}

inline void OneBodyDensityMatrices::normalizeBasis(ParticleSet& pset_target)
{
  int ngrid = std::max(200, input_.get_points());
  int ngtot = pow(ngrid, OHMMS_DIM);
  Real du   = input_.get_scale() / ngrid;
  Real dV   = volume_ / ngtot;
  Position rp;
  Vector<Value> bnorms;
  int gdims[OHMMS_DIM];
  gdims[0] = pow(ngrid, OHMMS_DIM - 1);
  for (int d = 1; d < OHMMS_DIM; ++d)
    gdims[d] = gdims[d - 1] / ngrid;
  bnorms.resize(basis_size_);
  for (int i = 0; i < basis_size_; ++i)
    bnorms[i] = 0.0;
  std::fill(basis_norms_.begin(), basis_norms_.end(), 1.0);
  for (int p = 0; p < ngtot; ++p)
  {
    int nrem = p;
    for (int d = 0; d < OHMMS_DIM - 1; ++d)
    {
      int ind = nrem / gdims[d];
      rp[d]   = ind * du + du / 2;
      nrem -= ind * gdims[d];
    }
    rp[OHMMS_DIM - 1] = nrem * du + du / 2;
    rp                = lattice_.toCart(rp) + rcorner_;
    updateBasis(rp, pset_target);
    for (int i = 0; i < basis_size_; ++i)
      bnorms[i] += qmcplusplus::conj(basis_values_[i]) * basis_values_[i] * dV;
  }
  for (int i = 0; i < basis_size_; ++i)
    basis_norms_[i] = 1.0 / std::sqrt(real(bnorms[i]));
}

void OneBodyDensityMatrices::registerOperatorEstimator(hdf_archive& file)
{
  std::vector<int> my_indexes(2, basis_size_);
  if constexpr (IsComplex_t<Value>::value)
  {
    my_indexes.push_back(2);
  }
  int nentries = std::accumulate(my_indexes.begin(), my_indexes.end(), 1);

  int spin_data_size = 0;
  if constexpr (IsComplex_t<Value>::value)
    spin_data_size = 2 * basis_size_ * basis_size_;
  else
    spin_data_size = basis_size_ * basis_size_;

  hdf_path hdf_name{my_name_};
  hdf_name /= "number_matrix";
  for (int s = 0; s < species_.size(); ++s)
  {
    h5desc_.emplace_back(hdf_name / species_.speciesName[s]);
    auto& oh = h5desc_.back();
    oh.set_dimensions(my_indexes, s * spin_data_size);
  }
}

} // namespace qmcplusplus
