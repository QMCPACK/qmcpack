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
#include <cmath>
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

  int nspecies        = species.size();
  if (! species_.hasAttribute("membersize"))
    throw UniformCommunicateError("OneBodyDensityMatrices::OneBodyDensityMatrices error: Species set does not have the required attribute 'membersize'");
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
  integrated_values.resize(basis_size);
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
  if(!input_.get_normalized())
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

void OneBodyDensityMatrices::warmupSampling(ParticleSet& pset_target, RandomGenerator_t& rng)
{
  if (sampling_ == Sampling::METROPOLIS)
  {
    if (!warmed_up)
    {
      rpcur = diffusion(std::sqrt(input_.get_timestep()), rng);
      rpcur += center_;
      if (input_.get_integrator() == Integrator::DENSITY)
        density_drift(rpcur, rhocur, dpcur, pset_target);
      else
        APP_ABORT("OneBodyDensityMatrices::warmup_sampling invalid integrator");
    }
    generateSamples(1.0, pset_target, rng, input_.get_warmup_samples());
    warmed_up = true;
  }
}


OneBodyDensityMatrices::FullPrecReal OneBodyDensityMatrices::evaluate_matrix(ParticleSet& P,
                                                                             MCPWalker& walker,
                                                                             ParticleSet& pset_target,
                                                                             TrialWaveFunction& psi_target,
                                                                             RandomGenerator_t& rng)
{
  warmupSampling(pset_target, rng);
  // get weight and single particle energy trace data
  Real weight;
  weight = walker.Weight * metric_;

  // compute sample positions (monte carlo or deterministic)
  generateSamples(weight, pset_target, rng);
  // compute basis and wavefunction ratio values in matrix form
  generate_sample_basis(Phi_MB, pset_target, psi_target);  // basis           : samples   x basis_size
  generate_sample_ratios(Psi_NM, pset_target, psi_target); // conj(Psi ratio) : particles x samples
  generate_particle_basis(P, Phi_NB, pset_target);         // conj(basis)     : particles x basis_size
  // perform integration via matrix products
  {
    ScopedTimer local_mp_timer(timers_.matrix_products_timer);
    for (int s = 0; s < species_.size(); ++s)
    {
      Matrix<Value>& Psi_nm     = Psi_NM[s];
      Matrix<Value>& Phi_Psi_nb = Phi_Psi_NB[s];
      Matrix<Value>& Phi_nb     = Phi_NB[s];
      diag_product(Psi_nm, sample_weights, Psi_nm);
      product(Psi_nm, Phi_MB, Phi_Psi_nb);       // ratio*basis : particles x basis_size
      product_AtB(Phi_nb, Phi_Psi_nb, N_BB[s]); // conj(basis)^T*ratio*basis : basis_size^2
    }
  }
  // accumulate data into collectables
  {
    ScopedTimer local_timer(timers_.accumulate_timer);
    const int basis_size2 = basis_size * basis_size;
    int ij                = nindex;
    for (int s = 0; s < species_.size(); ++s)
    {
      //int ij=nindex; // for testing
      const Matrix<Value>& NDM = N_BB[s];
      for (int n = 0; n < basis_size2; ++n)
      {
        Value val = NDM(n);
        P.Collectables[ij] += real(val);
        ij++;
#if defined(QMC_COMPLEX)
        P.Collectables[ij] += imag(val);
        ij++;
#endif
      }
    }
  }


#ifdef DMCHECK
  report();
  app_log() << "DM Check" << std::endl;
  evaluate_check(P);
  compare("  Phi_MB", Phi_MB, Phi_MBtmp);
  for (int s = 0; s < nspecies; ++s)
  {
    app_log() << "  species " << s << std::endl;
    compare("    E_N       ", *E_N[s], *E_Ntmp[s]);
    compare("    Phi_NB    ", *Phi_NB[s], *Phi_NBtmp[s]);
    compare("    Psi_NM    ", *Psi_NM[s], *Psi_NMtmp[s]);
    compare("    Phi_Psi_NB", *Phi_Psi_NB[s], *Phi_Psi_NBtmp[s], true);
    compare("    N_BB      ", *N_BB[s], *N_BBtmp[s], true);
  }
  app_log() << "end DM Check" << std::endl;
  APP_ABORT("DM Check");
#endif


  return 0.0;
}


OneBodyDensityMatrices::FullPrecReal OneBodyDensityMatrices::evaluate_check(ParticleSet& P)
{
#ifdef DMCHECK
  APP_ABORT("OneBodyDensityMatrices::evaluate_check  use of E_trace in this function needs to be replaces with "
            "get_energies() and E_samp");
  int n = 0;
  for (int s = 0; s < nspecies; ++s)
  {
    Matrix_t& Phi_mb     = Phi_MBtmp;
    Matrix_t& Psi_nm     = *Psi_NMtmp[s];
    Matrix_t& Phi_Psi_nb = *Phi_Psi_NBtmp[s];
    Matrix_t& Phi_nb     = *Phi_NBtmp[s];
    Vector_t& E_n        = *E_Ntmp[s];
    Matrix_t& N_bb       = *N_BBtmp[s];
    Matrix_t& E_bb       = *E_BBtmp[s];

    for (int ij = 0; ij < basis_size * basis_size; ++ij)
      N_bb(ij) = 0.0;
    for (int ij = 0; ij < basis_size * basis_size; ++ij)
      E_bb(ij) = 0.0;

    for (int ns = 0; ns < species_size[s]; ++ns, ++n)
    {
      std::fill(integrated_values.begin(), integrated_values.end(), 0.0);
      for (int m = 0; m < samples; ++m)
      {
        PosType& rsamp = rsamples[m];
        update_basis(rsamp, P);
        PosType dr = rsamp - P.R[n];
        P.makeMove(n, dr);
        Value_t ratio = sample_weights[m] * qmcplusplus::conj(Psi.calcRatio(P, n));
        P.rejectMove(n);
        for (int i = 0; i < basis_size; ++i)
        {
          integrated_values[i] += ratio * basis_values[i];
          Phi_mb(m, i) = basis_values[i];
        }
        Psi_nm(ns, m) = ratio;
      }
      update_basis(P.R[n]);
      for (int i = 0; i < basis_size; ++i)
        Phi_Psi_nb(ns, i) = integrated_values[i];
      for (int i = 0; i < basis_size; ++i)
        Phi_nb(ns, i) = qmcplusplus::conj(basis_values[i]);
      for (int i = 0; i < basis_size; ++i)
      {
        Value_t phi_i = qmcplusplus::conj(basis_values[i]);
        for (int j = 0; j < basis_size; ++j)
        {
          Value_t val = phi_i * integrated_values[j];
          N_bb(i, j) += val;
        }
      }
    }
  }
#endif
  return 0.0;
}


OneBodyDensityMatrices::FullPrecReal OneBodyDensityMatrices::evaluateLoop(ParticleSet& pset_target,
                                                                          MCPWalker& walker,
                                                                          TrialWaveFunction& psi_target,
                                                                          RandomGenerator_t& rng)
{
  const int basis_size2 = basis_size * basis_size;
  warmupSampling(pset_target, rng);
  Real weight;
  weight = walker.Weight * metric_;
  generateSamples(weight, pset_target, rng);
  int n = 0;
  for (int s = 0; s < species_.size(); ++s)
  {
    for (int ns = 0; ns < species_size[s]; ++ns, ++n)
    {
      integrate(pset_target, psi_target, n);
      update_basis(pset_target.R[n], pset_target);
      int ij = nindex + s * basis_size2;
      for (int i = 0; i < basis_size; ++i)
      {
        Value phi_i = qmcplusplus::conj(basis_values[i]);
        for (int j = 0; j < basis_size; ++j)
        {
          Value val = phi_i * integrated_values[j];
          pset_target.Collectables[ij] += real(val);
          ij++;
#if defined(QMC_COMPLEX)
          pset_target.Collectables[ij] += imag(val);
          ij++;
#endif
        }
      }
    }
  }
  return 0.0;
}

template<class RAN_GEN>
void OneBodyDensityMatrices::generateSamples(Real weight, ParticleSet& pset_target, RAN_GEN& rng, int steps)
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
    // trial pos
    Position rp;
    // trial drift
    Position dp;
    // drift sum
    Position ds;
    Real rhop, ratio, Pacc;              //trial density, dens ratio, acc prob
    Position diff = diffusion(sqt, rng); //get diffusion
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
{
  ScopedTimer scope_timer(timers_.eval_timer);
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    auto& pset_target = psets[iw];
    if (check_derivatives)
      test_derivatives(pset_target, rng);
    if (input_.get_evaluator() == Evaluator::LOOP)
      evaluateLoop(psets[iw], walkers[iw], wfns[iw], rng);
  }
  //   else if (input_.get_evaluator() == matrix)
  //     evaluate_matrix(P);
  //   else
  //     throw logic_error("OneBodyDensityMatrices::accumulate called with invalid evaluator, developer error");
  // }
  // return 0.0;
}


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


void OneBodyDensityMatrices::generate_sample_ratios(std::vector<Matrix<Value>>& Psi_nm,
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
    for (int s = 0; s < species_.size(); ++s)
    {
      Matrix<Value>& P_nm = Psi_nm[s];
      for (int n = 0; n < species_size[s]; ++n, ++p)
      {
        P_nm(n, m) = qmcplusplus::conj(psi_ratios[p]);
      }
    }
  }
}


void OneBodyDensityMatrices::generate_particle_basis(ParticleSet& P,
                                                     std::vector<Matrix<Value>>& Phi_nb,
                                                     ParticleSet& pset_target)
{
  ScopedTimer t(timers_.gen_particle_basis_timer);
  int p = 0;
  for (int s = 0; s < species_.size(); ++s)
  {
    int nb              = 0;
    Matrix<Value>& P_nb = Phi_nb[s];
    for (int n = 0; n < species_size[s]; ++n, ++p)
    {
      update_basis(P.R[p], pset_target);
      for (int b = 0; b < basis_size; ++b, ++nb)
        P_nb(nb) = qmcplusplus::conj(basis_values[b]);
    }
  }
}


inline void OneBodyDensityMatrices::integrate(ParticleSet& pset_target, TrialWaveFunction& psi_target, int n)
{
  std::fill(integrated_values.begin(), integrated_values.end(), 0.0);
  for (int s = 0; s < samples_; ++s)
  {
    Position& rsamp = rsamples[s];
    update_basis(rsamp, pset_target);
    pset_target.makeMove(n, rsamp - pset_target.R[n]);
    Value ratio = sample_weights[s] * qmcplusplus::conj(psi_target.calcRatio(pset_target, n));
    pset_target.rejectMove(n);
    for (int i = 0; i < basis_size; ++i)
      integrated_values[i] += ratio * basis_values[i];
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
  normalized = true;
}


void OneBodyDensityMatrices::test_derivatives(ParticleSet& pset_target, RandomGenerator_t& rng)
{
  app_log() << "OneBodyDensityMatrices::test_derivatives  checking drift" << std::endl;

  Position r, rtmp;

  Real delta = 1e-5;

  Real dens, densp, densm;
  Position drift, driftn, drifttmp;

  app_log() << "  warming up" << std::endl;
  warmupSampling(pset_target, rng);
  app_log() << "  generating samples" << std::endl;
  generateSamples(1.0, pset_target, rng);

  app_log() << "  testing derivatives at sample points" << std::endl;
  for (int s = 0; s < rsamples.size(); ++s)
  {
    r = rsamples[s];

    density_drift(r, dens, drift, pset_target);

    for (int d = 0; d < OHMMS_DIM; ++d)
    {
      rtmp = r;

      rtmp[d] = r[d] + delta;
      density_drift(rtmp, densp, drifttmp, pset_target);

      rtmp[d] = r[d] - delta;
      density_drift(rtmp, densm, drifttmp, pset_target);

      driftn[d] = (densp - densm) / (2 * delta);
    }
    driftn *= .5 * input_.get_timestep() / dens;

    app_log() << s << std::endl;
    app_log() << "  " << driftn << std::endl;
    app_log() << "  " << drift << std::endl;
  }

  APP_ABORT("OneBodyDensityMatrices::test_derivatives");
}

inline void OneBodyDensityMatrices::test_overlap(ParticleSet& pset_target)
{
  /** I feel like we should give the user what they want or abort.
   *  but is there an expectation that this test depends on the user set points?
   */
  int ngrid = std::max(50, input_.get_points());
  int ngtot = pow(ngrid, OHMMS_DIM);

  Position rp;
  Real du = input_.get_scale() / ngrid;
  Real dV = volume_ / ngtot;

  Position rmin = std::numeric_limits<Real>::max();
  Position rmax = -std::numeric_limits<Real>::max();
  int gdims[OHMMS_DIM];
  gdims[0] = pow(ngrid, OHMMS_DIM - 1);
  for (int d = 1; d < OHMMS_DIM; ++d)
    gdims[d] = gdims[d - 1] / ngrid;

  Array<Value, 2> omat;
  omat.resize(basis_size, basis_size);
  for (int i = 0; i < basis_size; ++i)
    for (int j = 0; j < basis_size; ++j)
      omat(i, j) = 0.0;

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
      for (int j = 0; j < basis_size; ++j)
        omat(i, j) += qmcplusplus::conj(basis_values[i]) * basis_values[j] * dV;
    for (int d = 0; d < OHMMS_DIM; ++d)
    {
      rmin[d] = std::min(rmin[d], rp[d]);
      rmax[d] = std::max(rmax[d], rp[d]);
    }
  }

  app_log() << "OneBodyDensityMatrices::test_overlap  checking overlap matrix" << std::endl;
  app_log() << "  rmin = " << rmin << std::endl;
  app_log() << "  rmax = " << rmax << std::endl;
  app_log() << "  overlap scale " << std::abs(omat(0, 0)) << std::endl;
  app_log() << "  overlap matrix:" << std::endl;
  for (int i = 0; i < basis_size; ++i)
  {
    app_log() << std::endl;
    for (int j = 0; j < basis_size; ++j)
      app_log() << std::abs(omat(i, j)) / std::abs(omat(0, 0)) << " ";
  }
  app_log() << std::endl;
  APP_ABORT("OneBodyDensityMatrices::test_overlap");
}


bool OneBodyDensityMatrices::match(Value e1, Value e2, Real tol)
{
  return std::abs(e1 - e2) < tol;
  //return std::abs(2*(e1-e2)/(e1+e2)) < tol;
}


bool OneBodyDensityMatrices::same(Vector<Value>& v1, Vector<Value>& v2, Real tol)
{
  if (v1.size() != v2.size())
    APP_ABORT("OneBodyDensityMatrices::same(vector)  vectors differ in size");
  bool sm = true;
  for (int i = 0; i < v1.size(); ++i)
    sm &= match(v1[i], v2[i]);
  return sm;
}

bool OneBodyDensityMatrices::same(Matrix<Value>& m1, Matrix<Value>& m2, Real tol)
{
  if (m1.rows() != m2.rows() || m1.cols() != m2.cols())
    APP_ABORT("OneBodyDensityMatrices::same(matrix)  matrices differ in size");
  bool sm = true;
  int n   = m1.rows() * m1.cols();
  for (int i = 0; i < n; ++i)
    sm &= match(m1(i), m2(i));
  return sm;
}

void OneBodyDensityMatrices::compare(const std::string& name,
                                     Vector<Value>& v1,
                                     Vector<Value>& v2,
                                     bool write,
                                     bool diff_only)
{
  bool sm            = same(v1, v2);
  std::string result = "differ";
  if (sm)
    result = "agree";
  app_log() << name << " " << result << std::endl;
  if (write && !sm)
    for (int i = 0; i < v1.size(); ++i)
      app_log() << "      " << i << " " << real(v1[i]) << " " << real(v2[i]) << " " << real(v1[i] / v2[i]) << " "
                << real(v2[i] / v1[i]) << std::endl;
}

void OneBodyDensityMatrices::compare(const std::string& name,
                                     Matrix<Value>& m1,
                                     Matrix<Value>& m2,
                                     bool write,
                                     bool diff_only)
{
  bool sm            = same(m1, m2);
  std::string result = "differ";
  if (sm)
    result = "agree";
  app_log() << name << " " << result << std::endl;
  if (write && !sm)
    for (int i = 0; i < m1.rows(); ++i)
      for (int j = 0; j < m1.cols(); ++j)
        if (!diff_only || !match(m1(i, j), m2(i, j)))
          app_log() << "      " << i << " " << j << " " << real(m1(i, j)) << " " << real(m2(i, j)) << " "
                    << real(m1(i, j) / m2(i, j)) << std::endl;
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
