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

OneBodyDensityMatrices::OneBodyDensityMatrices(const OneBodyDensityMatrices& obdm)
    : OperatorEstBase(obdm.data_locality_),
      input_(obdm.input_),
      lattice_(obdm.lattice_),
      basis_functions(obdm.basis_functions),
      timers_("OneBodyDensityMatrix")
{}

OneBodyDensityMatrices::OneBodyDensityMatrices(OneBodyDensityMatricesInput&& obdmi,
                                               const Lattice& lattice,
                                               const SpeciesSet& species)
    : OperatorEstBase(DataLocality::crowd), lattice_(lattice), species_(species), timers_("OneBodyDensityMatrix")
{}

OneBodyDensityMatrices::~OneBodyDensityMatrices() {}


OneBodyDensityMatrices* OneBodyDensityMatrices::clone() { return new OneBodyDensityMatrices(*this); }

void OneBodyDensityMatrices::startBlock(int steps) {}

void OneBodyDensityMatrices::setRandomGenerator(RandomGenerator_t* rng) { uniform_random = rng; }

inline void OneBodyDensityMatrices::generate_samples(Real weight, ParticleSet& pset_target, int steps)
{
  ScopedTimer local_timer(timers_.gen_samples_timer);
  RandomGenerator_t& rng = *uniform_random;
  bool save              = false;
  if (steps == 0)
  {
    save  = true;
    steps = samples_;
  }
  if (input_.get_integrator() == Integrator::UNIFORM_GRID)
    generate_uniform_grid(rng);
  else if (input_.get_integrator() == Integrator::UNIFORM)
    generate_uniform_samples(rng);
  else if (input_.get_integrator() == Integrator::DENSITY)
    generate_density_samples(save, steps, rng, pset_target);

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
  if (write_rstats && omp_get_thread_num() == 0)
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


inline void OneBodyDensityMatrices::generate_uniform_grid(RandomGenerator_t& rng)
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


inline void OneBodyDensityMatrices::generate_uniform_samples(RandomGenerator_t& rng)
{
  Position rp;
  for (int s = 0; s < samples_; ++s)
  {
    for (int d = 0; d < OHMMS_DIM; ++d)
      rp[d] = input_.get_scale() * rng();
    rsamples[s] = lattice_.toCart(rp) + rcorner_;
  }
}


inline void OneBodyDensityMatrices::generate_density_samples(bool save,
                                                             int steps,
                                                             RandomGenerator_t& rng,
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
    Position n, rp, dp, ds; //diffusion, trial pos/drift, drift sum
    Real rhop, ratio, Pacc; //trial density, dens ratio, acc prob
    diffusion(sqt, n);      //get diffusion
    if (input_.get_use_drift())
    {
      rp = r + n + d;                                                  //update trial position
      density_drift(rp, rhop, dp, pset_target);                        //get trial drift and density
      ratio = rhop / rho;                                              //density ratio
      ds    = dp + d;                                                  //drift sum
      Pacc  = ratio * std::exp(-ot * (dot(n, ds) + .5 * dot(ds, ds))); //acceptance probability
    }
    else
    {
      rp = r + n;                          //update trial position
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


inline void OneBodyDensityMatrices::diffusion(Real sqt, Position& diff)
{
  assignGaussRand(&diff[0], OHMMS_DIM, *uniform_random);
  diff *= sqt;
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


void OneBodyDensityMatrices::generate_sample_ratios(std::vector<Matrix<Value>*> Psi_nm,
                                                    ParticleSet& pset_target,
                                                    TrialWaveFunction& psi_target)
{
  ScopedTimer t(timers_.gen_sample_ratios_timer);
  for (int m = 0; m < samples_; ++m)
  {
    // get N ratios for the current sample point
    pset_target.makeVirtualMoves(rsamples[m]);
    psi_target.evaluateRatiosAlltoOne(pset_target, psi_ratios);

    // collect ratios into per-species matrices
    int p = 0;
    for (int s = 0; s < nspecies; ++s)
    {
      Matrix<Value>& P_nm = *Psi_nm[s];
      for (int n = 0; n < species_size[s]; ++n, ++p)
      {
        P_nm(n, m) = qmcplusplus::conj(psi_ratios[p]);
      }
    }
  }
}


void OneBodyDensityMatrices::generate_particle_basis(ParticleSet& P,
                                                     std::vector<Matrix<Value>*>& Phi_nb,
                                                     ParticleSet& pset_target)
{
  ScopedTimer t(timers_.gen_particle_basis_timer);
  int p = 0;
  for (int s = 0; s < nspecies; ++s)
  {
    int nb              = 0;
    Matrix<Value>& P_nb = *Phi_nb[s];
    for (int n = 0; n < species_size[s]; ++n, ++p)
    {
      update_basis(P.R[p], pset_target);
      for (int b = 0; b < basis_size; ++b, ++nb)
        P_nb(nb) = qmcplusplus::conj(basis_values[b]);
    }
  }
}

inline void OneBodyDensityMatrices::update_basis(const Position& r, ParticleSet& pset_target)
{
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

} // namespace qmcplusplus
