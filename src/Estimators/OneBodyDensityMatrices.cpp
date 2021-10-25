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
#include "QMCWaveFunctions/WaveFunctionFactory.h"


namespace qmcplusplus
{
using MatrixOperators::diag_product;
using MatrixOperators::product;
using MatrixOperators::product_AtB;

OneBodyDensityMatrices::OneBodyDensityMatrices(OneBodyDensityMatricesInput&& obdmi,
                                               const Lattice& lattice,
                                               const SpeciesSet& species,
                                               const WaveFunctionFactory& wf_factory)
    : OperatorEstBase(DataLocality::crowd),
      input_(obdmi),
      lattice_(lattice),
      species_(species),
      wf_factory_(wf_factory),
      timers_("OneBodyDensityMatrix")
{
  lattice_.reset();
  if (input_.get_center_defined())
    center_ = input_.get_center();
  else
    center_ = lattice_.Center;

  volume_   = lattice_.Volume * std::exp(OHMMS_DIM * std::log(input_.get_scale()));
  periodic_ = lattice_.SuperCellEnum != SUPERCELL_OPEN;
  rcorner_  = center_ - input_.get_scale() * lattice_.Center;

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
  rsamples.resize(samples_);

  // get the sposets that form the basis
  auto& sposets = input_.get_basis_sets();

  for (int i = 0; i < sposets.size(); ++i)
  {
    SPOSet* sposet = wf_factory.getSPOSet(sposets[i]);
    if (sposet == 0)
      throw UniformCommunicateError("OneBodyDensityMatrices::OneBodyDensityMatrices sposet " + sposets[i] +
                                    " does not exist");
    basis_functions.add(sposet->makeClone());
  }
  basis_size = basis_functions.size();

  if (basis_size < 1)
    throw UniformCommunicateError("OneBodyDensityMatrices::OneBodyDensityMatrices basis_size must be greater than one");

  int nspecies = species.size();
  if (!species_.hasAttribute("membersize"))
    throw UniformCommunicateError("OneBodyDensityMatrices::OneBodyDensityMatrices error: Species set does not have the "
                                  "required attribute 'membersize'");
  int isize = species.getAttribute("membersize");
  // We have the count per species at least a fundamental as the  total particles.
  // the sume of species membersize and total particles should be an invariant.
  int nparticles = 0;
  for (int s = 0; s < nspecies; ++s)
    nparticles += species(isize, s);
  for (int s = 0; s < nspecies; ++s)
    species_size.push_back(species(isize, s));
  for (int s = 0; s < nspecies; ++s)
    species_name.push_back(species.speciesName[s]);

  basis_values.resize(basis_size);
  basis_norms.resize(basis_size);

  Real bn_standard = 1.0;
  if (input_.get_volume_normalized())
    bn_standard = 1.0 / std::sqrt(volume_);
  for (int i = 0; i < basis_size; ++i)
    basis_norms[i] = bn_standard;

  rsamples.resize(samples_);
  sample_weights.resize(samples_);
  psi_ratios.resize(nparticles);

  if (input_.get_evaluator() == Evaluator::MATRIX)
  {
    Phi_MB.resize(samples_, basis_size);
    Phi_NB.reserve(nspecies);
    Psi_NM.reserve(nspecies);
    Phi_Psi_NB.reserve(nspecies);
    N_BB.reserve(nspecies);
    for (int s = 0; s < nspecies; ++s)
    {
      int specs_size = species_size[s];
      Phi_NB.emplace_back(specs_size, basis_size);
      Psi_NM.emplace_back(specs_size, samples_);
      Phi_Psi_NB.emplace_back(specs_size, basis_size);
      N_BB.emplace_back(basis_size, basis_size);
    }
  }

  if (sampling_ == Sampling::METROPOLIS)
  {
    basis_gradients.resize(basis_size);
    basis_laplacians.resize(basis_size);
  }

  // so if the input is not normalized, normalize it.
  // with respect to what?
  if (!input_.get_normalized())
  {
    auto& pset_target = wf_factory_.getTargetParticleSet();
    normalize(pset_target);
  }
}

OneBodyDensityMatrices::OneBodyDensityMatrices(const OneBodyDensityMatrices& obdm)
    : OneBodyDensityMatrices(OneBodyDensityMatricesInput(obdm.input_), obdm.lattice_, obdm.species_, obdm.wf_factory_)
{}

OneBodyDensityMatrices::~OneBodyDensityMatrices() {}

OneBodyDensityMatrices* OneBodyDensityMatrices::clone() { return new OneBodyDensityMatrices(*this); }

void OneBodyDensityMatrices::startBlock(int steps) {}

template<class RNG_GEN>
void OneBodyDensityMatrices::generateSamples(Real weight, ParticleSet& pset_target, RNG_GEN& rng, int steps)
{
  ScopedTimer local_timer(timers_.gen_samples_timer);

  // \todo make this an ifndef NDEBUG section which is much clearer.
  bool save = false;
  if (steps == 0)
  {
    save  = true;
    steps = samples_;
  }

  switch (input_.get_integrator())
  {
  case Integrator::UNIFORM_GRID:
    generate_uniform_grid(rng);
    break;
  case Integrator::UNIFORM:
    generate_uniform_samples(rng);
    break;
  case Integrator::DENSITY:
    generate_density_samples(save, steps, rng, pset_target);
    break;
  }

  if (save)
  {
    if (sampling_ == Sampling::METROPOLIS)
    {
      sample_weights *= weight;
    }
    else
    {
      std::fill(sample_weights.begin(), sample_weights.end(), weight);
    }
  }

  // temporary check
  if (input_.get_rstats() && omp_get_thread_num() == 0)
  {
    Position rmin  = std::numeric_limits<Real>::max();
    Position rmax  = -std::numeric_limits<Real>::max();
    Position rmean = 0.0;
    Position rstd  = 0.0;
    for (int s = 0; s < rsamples.size(); ++s)
      for (int d = 0; d < OHMMS_DIM; ++d)
      {
        Real rd = rsamples[s][d];
        rmin[d] = std::min(rmin[d], rd);
        rmax[d] = std::max(rmax[d], rd);
        rmean[d] += rd;
        rstd[d] += rd * rd;
      }
    rmean /= rsamples.size();
    rstd /= rsamples.size();
    for (int d = 0; d < OHMMS_DIM; ++d)
      rstd[d] = std::sqrt(rstd[d] - rmean[d] * rmean[d]);
    app_log() << "\nrsamples properties:" << std::endl;
    app_log() << "  rmin  = " << rmin << std::endl;
    app_log() << "  rmax  = " << rmax << std::endl;
    app_log() << "  rmean = " << rmean << std::endl;
    app_log() << "  rstd  = " << rstd << std::endl;
  }
}

template<typename RNG_GEN>
inline void OneBodyDensityMatrices::generate_uniform_grid(RNG_GEN& rng)
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
    rsamples[s]       = lattice_.toCart(rp) + rcorner_;
  }
}

template<typename RAN_GEN>
inline void OneBodyDensityMatrices::generate_uniform_samples(RAN_GEN& rng)
{
  Position rp;
  for (int s = 0; s < samples_; ++s)
  {
    for (int d = 0; d < OHMMS_DIM; ++d)
      rp[d] = input_.get_scale() * rng();
    rsamples[s] = lattice_.toCart(rp) + rcorner_;
  }
}

template<typename RAN_GEN>
inline void OneBodyDensityMatrices::generate_density_samples(bool save,
                                                             int steps,
                                                             RAN_GEN& rng,
                                                             ParticleSet& pset_target)
{
  const auto timestep = input_.get_timestep();
  Real sqt            = std::sqrt(timestep);
  Real ot             = 1.0 / timestep;
  Position r          = rpcur;  //current position
  Position d          = dpcur;  //current drift
  Real rho            = rhocur; //current density
  for (int s = 0; s < steps; ++s)
  {
    nmoves++;
    Position rp;                         // trial pos
    Position dp;                         // trial drift
    Position ds;                         // drift sum
    Real rhop;                           // trial density
    Real ratio;                          // dens ratio
    Real Pacc;                           // acc prob
    Position diff = diffusion(sqt, rng); // get diffusion
    if (input_.get_use_drift())
    {
      rp = r + diff + d;                                                  //update trial position
      density_drift(rp, rhop, dp, pset_target);                           //get trial drift and density
      ratio = rhop / rho;                                                 //density ratio
      ds    = dp + d;                                                     //drift sum
      Pacc  = ratio * std::exp(-ot * (dot(diff, ds) + .5 * dot(ds, ds))); //acceptance probability
    }
    else
    {
      rp = r + diff;                       //update trial position
      density_only(rp, rhop, pset_target); //get trial density
      ratio = rhop / rho;                  //density ratio
      Pacc  = ratio;                       //acceptance probability
    }
    if (rng() < Pacc)
    { //accept move
      r   = rp;
      d   = dp;
      rho = rhop;
      naccepted++;
    }
    if (save)
    {
      rsamples[s]       = r;
      sample_weights[s] = 1.0 / rho;
    }
  }
  acceptance_ratio = Real(naccepted) / nmoves;

  if (write_acceptance_ratio && omp_get_thread_num() == 0)
    app_log() << "dm1b  acceptance_ratio = " << acceptance_ratio << std::endl;

  rpcur  = r;
  dpcur  = d;
  rhocur = rho;
}

template<typename RAN_GEN>
OneBodyDensityMatrices::Position OneBodyDensityMatrices::diffusion(const Real sqt, RAN_GEN& rng)
{
  Position diff;
  assignGaussRand(&diff[0], OHMMS_DIM, rng);
  diff *= sqt;
  return diff;
}


inline void OneBodyDensityMatrices::density_only(const Position& r, Real& dens, ParticleSet& pset_target)
{
  update_basis(r, pset_target);
  dens = 0.0;
  for (int i = 0; i < basis_size; ++i)
  {
    Value b = basis_values[i];
    dens += std::abs(qmcplusplus::conj(b) * b);
  }
  dens /= basis_size;
}


void OneBodyDensityMatrices::density_drift(const Position& r, Real& dens, Position& drift, ParticleSet& pset_target)
{
  update_basis_d012(r, pset_target);
  dens  = 0.0;
  drift = 0.0;
  for (int i = 0; i < basis_size; ++i)
  {
    const Grad& bg = basis_gradients[i];
    Value b        = basis_values[i];
    Value bc       = qmcplusplus::conj(b);
    dens += std::abs(bc * b);
    for (int d = 0; d < OHMMS_DIM; ++d)
      drift[d] += std::real(bc * bg[d]);
  }
  drift *= input_.get_timestep() / dens;
  dens /= basis_size;
}

void OneBodyDensityMatrices::accumulate(const RefVector<MCPWalker>& walkers,
                                        const RefVector<ParticleSet>& psets,
                                        const RefVector<TrialWaveFunction>& wfns,
                                        RandomGenerator_t& rng)
{}

void OneBodyDensityMatrices::generate_sample_basis(Matrix<Value>& Phi_mb,
                                                   ParticleSet& pset_target,
                                                   TrialWaveFunction& psi_target)
{
  ScopedTimer local_timer(timers_.gen_sample_basis_timer);
  int mb = 0;
  for (int m = 0; m < samples_; ++m)
  {
    update_basis(rsamples[m], pset_target);
    for (int b = 0; b < basis_size; ++b, ++mb)
      Phi_mb(mb) = basis_values[b];
  }
}

inline void OneBodyDensityMatrices::update_basis(const Position& r, ParticleSet& pset_target)
{
  // This is ridiculous in the case of splines, still necessary for hybrid/LCAO
  pset_target.makeMove(0, r - pset_target.R[0]);
  basis_functions.evaluateValue(pset_target, 0, basis_values);
  pset_target.rejectMove(0);
  for (int i = 0; i < basis_size; ++i)
    basis_values[i] *= basis_norms[i];
}


inline void OneBodyDensityMatrices::update_basis_d012(const Position& r, ParticleSet& pset_target)
{
  pset_target.makeMove(0, r - pset_target.R[0]);
  basis_functions.evaluateVGL(pset_target, 0, basis_values, basis_gradients, basis_laplacians);
  pset_target.rejectMove(0);
  for (int i = 0; i < basis_size; ++i)
    basis_values[i] *= basis_norms[i];
  for (int i = 0; i < basis_size; ++i)
    basis_gradients[i] *= basis_norms[i];
  for (int i = 0; i < basis_size; ++i)
    basis_laplacians[i] *= basis_norms[i];
}

inline void OneBodyDensityMatrices::normalize(Real invToWgt) {}

inline void OneBodyDensityMatrices::normalize(ParticleSet& pset_target)
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
  bnorms.resize(basis_size);
  for (int i = 0; i < basis_size; ++i)
    bnorms[i] = 0.0;
  std::fill(basis_norms.begin(), basis_norms.end(), 1.0);
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
    update_basis(rp, pset_target);
    for (int i = 0; i < basis_size; ++i)
      bnorms[i] += qmcplusplus::conj(basis_values[i]) * basis_values[i] * dV;
  }
  for (int i = 0; i < basis_size; ++i)
    basis_norms[i] = 1.0 / std::sqrt(real(bnorms[i]));
}


template void OneBodyDensityMatrices::generateSamples<RandomGenerator_t>(Real weight,
                                                                         ParticleSet& pset_target,
                                                                         RandomGenerator_t& rng,
                                                                         int steps);
template void OneBodyDensityMatrices::generateSamples<StdRandom<double>>(Real weight,
                                                                         ParticleSet& pset_target,
                                                                         StdRandom<double>& rng,
                                                                         int steps);
} // namespace qmcplusplus
