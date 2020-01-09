//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef EWALD_TOOLS_H
#define EWALD_TOOLS_H

#include "QMCHamiltonians/EwaldRef.h"

namespace qmcplusplus
{
namespace ewaldtools
{
  using ewaldref::int_t;
  using ewaldref::real_t;
  using ewaldref::IntVec;
  using ewaldref::RealVec;
  using ewaldref::RealMat;
  using ewaldref::ChargeArray;
  using ewaldref::PosArray;
  using ewaldref::getKappaMadelung;
  using ewaldref::getKappaEwald;
  using ewaldref::RspaceMadelungTerm;
  using ewaldref::KspaceMadelungTerm;
  using ewaldref::RspaceEwaldTerm;
  using ewaldref::KspaceEwaldTerm;


/// Similar to ewaldref::gridSum but for adaptive anisotropic grids
template<typename T>
real_t anisotropicGridSum(T& function, IntVec& nmax, bool zero = true, real_t tol = 1e-11)
{
  real_t dv  = std::numeric_limits<real_t>::max();
  real_t dva = std::numeric_limits<real_t>::max();
  real_t v   = 0.0;
  int_t im   = 0;
  int_t jm   = 0;
  int_t km   = 0;
  IntVec iv;
  // Add the value at the origin, if requested.
  if (zero)
  {
    iv = 0;
    v += function(iv);
  }
  // Sum over cubic surface shells until the tolerance is reached.
  tol /= 3; // Tolerance is to be applied to each pair of surface planes separately
  int_t nx = 0;
  int_t ny = 0;
  int_t nz = 0;
  bool last_surface = false;
  while(!last_surface)
  {
    // Add surface planes until none contribute.
    // Directions are turned off selectively once they no longer contribute.
    bool x_active = true;
    bool y_active = true;
    bool z_active = true;
    int_t nx_added = 0;
    int_t ny_added = 0;
    int_t nz_added = 0;
    while (x_active || y_active || z_active)
    {
      dv = 0.0;
      // Sum over new surface planes perpendicular to the x direction.
      if(x_active)
      {
        dva = 0.0;
        im += 1;
        for (int_t i : {-im, im})
          for (int_t j = -jm; j < jm + 1; ++j)
            for (int_t k = -km; k < km + 1; ++k)
            {
              iv[0]    = i;
              iv[1]    = j;
              iv[2]    = k;
              real_t f = function(iv);
              dv += f;
              dva += std::abs(f);
            }
        nx_added++;
        x_active = dva > tol; // Stop sum in x-direction when not contributing.
      }
      // Sum over new surface planes perpendicular to the y direction.
      if(y_active)
      {
        dva = 0.0;
        jm += 1;
        for (int_t j : {-jm, jm})
          for (int_t k = -km; k < km + 1; ++k)
            for (int_t i = -im; i < im + 1; ++i)
            {
              iv[0]    = i;
              iv[1]    = j;
              iv[2]    = k;
              real_t f = function(iv);
              dv += f;
              dva += std::abs(f);
            }
        ny_added++;
        y_active = dva > tol; // Stop sum in y-direction when not contributing.
      }
      // Sum over new surface planes perpendicular to the z direction.
      if(z_active)
      {
        dva = 0.0;
        km += 1;
        for (int_t k : {-km, km})
          for (int_t i = -im; i < im + 1; ++i)
            for (int_t j = -jm; j < jm + 1; ++j)
            {
              iv[0]    = i;
              iv[1]    = j;
              iv[2]    = k;
              real_t f = function(iv);
              dv += f;
              dva += std::abs(f);
            }
        nz_added++;
        z_active = dva > tol; // Stop sum in z-direction when not contributing.
      }
      v += dv;
    }
    nx += nx_added;
    ny += ny_added;
    nz += nz_added;
    last_surface = nx_added==1 && ny_added==1 && nz_added==1;
  }

  // Return dimensions of the converged anisotropic grid
  nmax[0] = im;
  nmax[1] = jm;
  nmax[2] = km;

  return v;
}


/// Similar to ewaldref::gridSum, but for non-adaptive, fixed size grids
template<typename T>
real_t fixedGridSum(T& function, IntVec nmax, bool zero=true)
{
  real_t v = 0.0;
  IntVec iv;
  if(zero)
  {
    // Directly sum over the full space
    for(int_t i = -nmax[0]; i < nmax[0] + 1; ++i)
      for(int_t j = -nmax[1]; j < nmax[1] + 1; ++j)
        for(int_t k = -nmax[2]; k < nmax[2] + 1; ++k)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          v += function(iv);
        }
  }
  else
  {
    // Sum half spaces around zero in z
    for(int_t k = 1; k < nmax[2] + 1; ++k)
      for(int_t i = -nmax[0]; i < nmax[0] + 1; ++i)
        for(int_t j = -nmax[1]; j < nmax[1] + 1; ++j)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          v += function(iv);
          iv[2] = -k;
          v += function(iv);
        }
    // Sum half planes around zero in y
    iv[2] = 0;
    for(int_t j = 1; j < nmax[1] + 1; ++j)
      for(int_t i = -nmax[0]; i < nmax[0] + 1; ++i)
        {
          iv[0] = i;
          iv[1] = j;
          v += function(iv);
          iv[1] = -j;
          v += function(iv);
        }
    // Sum half lines around zero in x
    iv[1] = 0;
    for(int_t i = 1; i < nmax[0] + 1; ++i)
    {
      iv[0] = i;
      v += function(iv);
      iv[0] = -i;
      v += function(iv);
    }
  }
  return v;
}



/// Functor for short range part of e-e pair potential
class RspaceEwaldPairPotential
{
private:
  /// The constant 1/(\sqrt{2}\kappa) in Drummond 2008 formula 6
  const real_t rexponent;

public:
  RspaceEwaldPairPotential(real_t rexponent_in) : rexponent(rexponent_in) {}

  real_t operator()(const RealVec& r) const
  {
    real_t R  = std::sqrt(dot(r, r));
    return std::erfc(rexponent * R) / R;
  }
};


/// Functor for Fourier transform of long range part of e-e pair potential
class KspaceEwaldPairPotential
{
private:
  /// The constant -\kappa^2/2 in Drummond 2008 formula 6
  const real_t kexponent;
  /// The constant 4\pi/\Omega in Drummond 2008 formula 6
  const real_t kprefactor;

public:
  KspaceEwaldPairPotential(real_t kexponent_in, real_t kprefactor_in)
    : kexponent(kexponent_in), kprefactor(kprefactor_in)
  {}

  real_t operator()(const RealVec& k) const
  {
    real_t k2 = dot(k,k);
    return kprefactor * std::exp(kexponent * k2) / k2;
  }
};



template<typename M>
RealVec wignerPoint(const M& A, size_t nc=5)
{
  RealVec rwig = 0.0;
  real_t rmin = 1e90;
  for(int_t k=-nc; k<nc+1; k++)
    for(int_t j=-nc; j<nc+1; j++)
      for(int_t i=-nc; i<nc+1; i++)
        if(i!=0 || j!=0 || k!=0)
        {
          IntVec n(i,j,k);
          RealVec rv = dot(n,A);
          real_t rm = 0.5*std::sqrt(dot(rv,rv));
          if(rm<rmin)
          {
            rwig = rv;
            rmin = rm;
          }
        }
  return rwig;
}


template<typename M>
real_t wignerRadius(const M& A, size_t nc=5)
{
  const auto& rwig = wignerPoint(A,nc);
  return std::sqrt(dot(rwig,rwig));
}


template<typename V>
V ivmax(const V& v1, const V& v2)
{
  V vm;
  for(size_t d: {0,1,2})
    vm[d] = std::max(v1[d],v2[d]);
  return vm;
}



/// General class to perform Ewald sums
class AnisotropicEwald
{
public:
  IntVec nrmax; // Max aniso r-space grid dims found across all adaptive evals
  IntVec nkmax; // Max aniso k-space grid dims found across all adaptive evals

private:
  RealMat A;      // Real-space cell axes
  RealMat B;      // K-space cell axes
  ChargeArray Q;  // List of charges for all particles
  real_t tol;     // Tolerance in Ha for total energy adaptive sums

  size_t N;       // Number of particles
  real_t volume;  // Volume of real-space cell, \Omega in Drummond 2008

  // Constants used in Madelung sums (see Drummond 2008 formula 7)
  real_t kappa_madelung;      // kappa for Madelung sums
  real_t rexponent_madelung;  // The constant \kappa
  real_t kexponent_madelung;  // The constant 1/(4\kappa^2)
  real_t kprefactor_madelung; // The constant 4\pi/\Omega

  // Constants used in Ewald sums (see Drummond 2008 formula 6)
  real_t kappa_ewald;         // kappa for Ewald sums
  real_t rexponent_ewald;     // The constant 1/(\sqrt{2}\kappa)
  real_t kexponent_ewald;     // The constant -\kappa^2/2
  real_t kprefactor_ewald;    // The constant 4\pi/\Omega
  real_t kconstant_ewald;     // The constant -2*\pi*\kappa^2/\Omega

  // Internally stored Madelung values (computed upon instantiation)
  real_t madelung_constant;   // v_M in Drummond 2008 formula 5
  real_t madelung_energy;     // Total Madelung energy for all particles

  // Maximum anisotropic grid found across all adaptive evaluations
  IntVec nmax_anisotropic;

  // Internally stored values for optimized computation
  real_t ewald_constant_lr_energy; // Constant portion of LR (k-space) total energy
  std::vector<RealVec> kpoints;    // Converged k-point grid
  std::vector<real_t> vk;          // Fourier transform of Ewald pair potential
  real_t vk_sum;                   // Constant sum term in rho_k breakup
  std::vector<real_t> rhok_r;      // Charge structure factor
  std::vector<real_t> rhok_i;      // Charge structure factor

public:

  /// Empty constructor
  AnisotropicEwald() { }

  /// State-initializing constructor
  AnisotropicEwald(const RealMat& A_in, const ChargeArray& Q_in, real_t tol_in = 1e-10,real_t kappa_in=-1.0)
  {
    initialize(A_in,Q_in,tol_in,kappa_in);
  }

  /// Initialize constant data, including Madelung sums
  void initialize(const RealMat& A_in, const ChargeArray& Q_in, real_t tol_in = 1e-10,real_t kappa_in=-1.0)
  {
    // Zero out observed grid dimensions
    nrmax = 0;
    nkmax = 0;

    // Store input data
    A   = A_in;
    Q   = Q_in;
    tol = tol_in;

    // Setup cell data
    N      = Q.size();                         // Number of particles
    volume = std::abs(det(A));                 // Volume
    B      = 2 * M_PI * transpose(inverse(A)); // Reciprocal lattice vectors

    // Set Ewald sum constants
    // Madelung constants
    kappa_madelung      = getKappaMadelung(volume);
    rexponent_madelung  = kappa_madelung;
    kexponent_madelung  = -1. / (4 * std::pow(kappa_madelung, 2));
    kprefactor_madelung = 4 * M_PI / volume;
    // Ewald constants
    setKappaEwald(kappa_in);

    // Calculate and store Madelung energy
    madelung_energy = madelungEnergy();
  }


  void setKappaEwald(real_t kappa=-1.0)
  {
    if(kappa<0.0)
      kappa_ewald = getKappaEwald(volume);
    else
      kappa_ewald = kappa;
    rexponent_ewald  = 1. / (std::sqrt(2.) * kappa_ewald);
    kexponent_ewald  = -std::pow(kappa_ewald, 2) / 2;
    kprefactor_ewald = 4 * M_PI / volume;
    kconstant_ewald  = -2 * M_PI * std::pow(kappa_ewald, 2) / volume;

    // Calculate and store constant energy term for LR interaction
    real_t qsum  = 0.0;
    real_t qsum2 = 0.0;
    for(auto q: Q)
    {
      qsum  += q;
      qsum2 += q*q;
    }
    ewald_constant_lr_energy = 0.5*(qsum*qsum - qsum2)*kconstant_ewald;
  }


  /// Calculate the Madelung constant
  real_t madelungConstant() const
  {
    // Find maximum self-interaction charge product
    real_t QQmax = 0.0;
    for (auto q: Q)
      QQmax = std::max(q*q, QQmax);

    // Set Madelung tolerance
    real_t tol_madelung = tol * 2. / QQmax;

    // Create real-/k-space functors for terms within the Madelung sums
    RspaceMadelungTerm rfunc(A, rexponent_madelung);
    KspaceMadelungTerm kfunc(B, kexponent_madelung, kprefactor_madelung);
    IntVec nmax;

    // Compute the constant term
    real_t cval = -M_PI / (std::pow(kappa_madelung, 2) * volume) - 2 * kappa_madelung / std::sqrt(M_PI);
    // Compute the real-space sum (excludes zero)
    real_t rval = anisotropicGridSum(rfunc, nmax, false, tol_madelung);
    // Compute the k-space sum (excludes zero)
    real_t kval = anisotropicGridSum(kfunc, nmax, false, tol_madelung);

    // Sum all contributions to get the Madelung constant
    real_t vm = cval + rval + kval;

    return vm;
  }


  /// Calculate and store the Madelung energy
  real_t madelungEnergy()
  {
    // Madelung energy for a single particle with charge 1
    madelung_constant = madelungConstant();
    
    // Sum the Madelung self interaction for each particle
    real_t me = 0.0;
    for (size_t i = 0; i < N; ++i)
      me += Q[i] * Q[i] * madelung_constant / 2;

    return me;
  }



  /// Pair interaction computed via sum over adaptive anisotropic grid
  real_t ewaldPairPotential(const RealVec& r, real_t tol_ewald)
  {
    // Create real-/k-space functors for terms within the Ewald pair potential sums
    RspaceEwaldTerm rfunc(r, A, rexponent_ewald);
    KspaceEwaldTerm kfunc(r, B, kexponent_ewald, kprefactor_ewald);
    IntVec nmax;

    // Compute the constant term
    real_t cval = kconstant_ewald;
    // Compute the real-space sum (includes zero)
    real_t rval = anisotropicGridSum(rfunc, nmax, true, tol_ewald);
    nrmax = ivmax(nmax,nrmax);
    // Compute the k-space sum (excludes zero)
    real_t kval = anisotropicGridSum(kfunc, nmax, false, tol_ewald);
    nkmax = ivmax(nmax,nkmax);

    // Sum all contributions to get the Ewald pair interaction
    real_t epp = cval + rval + kval;

    return epp;
  }


  /// Total energy for all pairs using adaptive anisotropic grids
  template<typename PA>
  real_t ewaldEnergy(const PA& R)
  {
    nrmax = 0;
    nkmax = 0;
    real_t ee = madelung_energy;
    size_t Npairs = (N*(N-1))/2;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
        ee += Q[i]*Q[j]*ewaldPairPotential(R[i]-R[j],tol/(Q[i]*Q[j])/Npairs);
    return ee;
  }


  /// Total energy for all pairs using fixed anisotropic grid
  template<typename PA>
  real_t fixedGridEwaldEnergy(const PA& R, const IntVec& nmax)
  {
    real_t ee = madelung_energy;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
      {
        RealVec r = R[i]-R[j];

        // Create real-/k-space fuctors for terms within the Ewald pair potential sums
        RspaceEwaldTerm rfunc(r, A, rexponent_ewald);
        KspaceEwaldTerm kfunc(r, B, kexponent_ewald, kprefactor_ewald);
    
        // Compute the constant term
        real_t cval = kconstant_ewald;
        // Compute the real-space sum (includes zero)
        real_t rval = fixedGridSum(rfunc, nmax, true);
        // Compute the k-space sum (excludes zero)
        real_t kval = fixedGridSum(kfunc, nmax, false);

        // Sum all contributions to get the Ewald pair interaction
        real_t epp = cval + rval + kval;

        ee += Q[i]*Q[j]*epp;
      }
    return ee;
  }


  /// Total energy for all pairs using fixed maximum anisotropic grid
  template<typename PA>
  real_t fixedGridEwaldEnergy(const PA& R)
  {
    return fixedGridEwaldEnergy(R, nmax_anisotropic);
  }


  /// Constant part of total energy
  real_t ewaldEnergyConst()
  {
    return madelung_energy;
  }


  /// SR (real-space) part of total energy computed adaptively 
  template<typename PA>
  real_t ewaldEnergySR(const PA& R)
  {
    nrmax = 0;
    real_t ee = 0.0;
    size_t Npairs = (N*(N-1))/2;
    IntVec nmax;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
      {
        real_t tol_ewald = tol/(Q[i]*Q[j])/Npairs;
        RspaceEwaldTerm rfunc(R[i]-R[j], A, rexponent_ewald);
        real_t rval = anisotropicGridSum(rfunc, nmax, true, tol_ewald);
        nrmax = ivmax(nmax,nrmax);
        ee += Q[i]*Q[j]*rval;
      }
    return ee;
  }


  /// LR (k-space) part of total energy computed adaptively 
  template<typename PA>
  real_t ewaldEnergyLR(const PA& R)
  {
    nkmax = 0;
    real_t ee = ewald_constant_lr_energy;
    size_t Npairs = (N*(N-1))/2;
    IntVec nmax;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
      {
        real_t tol_ewald = tol/(Q[i]*Q[j])/Npairs;
        KspaceEwaldTerm kfunc(R[i]-R[j], B, kexponent_ewald, kprefactor_ewald);
        real_t kval = anisotropicGridSum(kfunc, nmax, false, tol_ewald);
        nkmax = ivmax(nmax,nkmax);
        ee += Q[i]*Q[j]*kval;
      }
    return ee;
  }


  
  /// Setup data structures for optimized computation on a fixed anisotropic k-grid
  void setupOpt(IntVec nmax)
  {
    // Calculate charge product sum
    real_t QQsum = 0.0;
    for(auto q: Q)
      QQsum += q*q;

    // Store k-points and Fourier transform of LR part of Ewald pair potential
    IntVec grid = 2*nmax+1;
    size_t nkpoints = grid[0]*grid[1]*grid[2]-1;
    kpoints.resize(nkpoints);
    vk.resize(nkpoints);
    vk_sum = 0.0;
    KspaceEwaldPairPotential vkf(kexponent_ewald,kprefactor_ewald);
    size_t n=0;
    for(int_t i = -nmax[0]; i < nmax[0] + 1; ++i)
      for(int_t j = -nmax[1]; j < nmax[1] + 1; ++j)
        for(int_t k = -nmax[2]; k < nmax[2] + 1; ++k)
        {
          if(i==0 && j==0 && k==0)
            continue;
          IntVec iv(i,j,k);
          RealVec kp = dot(iv,B);
          kpoints[n] = kp;
          real_t vk_val = vkf(kp);
          vk[n]      = vk_val;
          vk_sum    += 0.5*QQsum*vk_val;
          n++;
        }

    rhok_r.resize(nkpoints);
    rhok_i.resize(nkpoints);
  }


  // Computes error bound for excluding all images of SR sum for single e-e interaction
  real_t SROptErrorBound(RealVec wigner_point,real_t kappa,real_t errtol=1e-10)
  {
    real_t rexponent = 1./(std::sqrt(2.)*kappa);
    RspaceEwaldTerm vsr(wigner_point,A,rexponent);
    IntVec nmax;
    real_t error_bound =0.0;//= ansitropicGridSum(vsr,nmax,false,errtol);
    return error_bound;
  }


  // Find kappa to use in optmized computation based on error tolerance of SR part
  real_t findOptKappa(real_t errtol)
  {
    real_t qsum  = 0.0;
    real_t qsum2 = 0.0;
    for(auto q: Q)
    {
      qsum  += std::abs(q);
      qsum2 += q*q;
    }
    real_t QQbound = 0.5*(qsum*qsum - qsum2);
    RealVec wigner_point = wignerPoint(A);
    real_t wigner_radius = std::sqrt(dot(wigner_point,wigner_point));
    real_t err_bound_tol = errtol/QQbound/1e2;
    real_t error = 1e90;
    real_t kappa = wigner_radius;
    int_t n = 0;
    while(error>errtol && kappa>0)
    {
      kappa = (1.0-0.1*n)*wigner_radius;
      error = QQbound*SROptErrorBound(wigner_point,kappa);
    }
  }


  template<typename PA, typename DT>
  real_t ewaldEnergyOpt(const PA& R, const DT& dt)
  {
    real_t Vc  = ewaldEnergyConst();
    real_t Vsr = ewaldEnergySROpt(dt);
    real_t Vlr = ewaldEnergyLROpt(R);
    return Vc + Vsr + Vlr;
  }


  /// SR (r-space) part of total energy computed optimally
  template<typename DT>
  real_t ewaldEnergySROpt(const DT& dt)
  {
    const auto& dr = dt.getDisplacements();
    real_t ee = 0.0;
    RspaceEwaldPairPotential vsr(rexponent_ewald);
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
        ee += Q[i]*Q[j]*vsr(dr[i][j]);
    return ee;
  }

  
  /// LR (k-space) part of total energy computed optimally
  template<typename PA>
  real_t ewaldEnergyLROpt(const PA& R)
  {
    real_t ee = ewald_constant_lr_energy;

    real_t ve = 0.0;
//#pragma omp parallel for reduction(+ : ve)
    for(size_t n=0; n<kpoints.size(); n++)
    {
      const auto& k = kpoints[n];
      real_t coskr = 0.0;
      real_t sinkr = 0.0;
//#pragma omp simd
      for(size_t i=0; i<N; i++)
      {
        real_t c,s;
        real_t q = Q[i];
        sincos(dot(k,R[i]),&s,&c);
        coskr += q*c;
        sinkr += q*s;
      }
      ve += (coskr*coskr + sinkr*sinkr)*vk[n];
    }

    ee += 0.5*ve - vk_sum;

    return ee;
  }

  
  /// LR (k-space) part of total energy computed optimally
  template<typename PA>
  real_t ewaldEnergyLROpt_qmcpack(const PA& R)
  {
    real_t ee = ewald_constant_lr_energy;

    std::fill(rhok_r.begin(),rhok_r.end(),0.0);
    std::fill(rhok_i.begin(),rhok_i.end(),0.0);
    for(size_t i=0; i<N; i++)
    {
      real_t q = Q[i];
      const auto& r = R[i];
      for(size_t n=0; n<kpoints.size(); n++)
      {
        real_t c,s;
        sincos(dot(kpoints[n],r),&s,&c);
        rhok_r[n] += q*c;
        rhok_i[n] += q*s;
      }
    }

    real_t ve = 0.0;
    for(size_t n=0; n<kpoints.size(); n++)
    {
      ve += (rhok_r[n]*rhok_r[n]+rhok_i[n]*rhok_i[n])*vk[n];
    }

    ee += 0.5*ve - vk_sum;

    return ee;
  }


  /// SR (r-space) part of total energy computed optimally
  template<typename DT>
  real_t ewaldEnergySROpt_orig(const DT& dt)
  {
    const auto& dr = dt.getDisplacements();
    real_t ee = 0.0;
    RspaceEwaldPairPotential vsr(rexponent_ewald);
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
        ee += Q[i]*Q[j]*vsr(dr[i][j]);
    return ee;
  }

  
  /// LR (k-space) part of total energy computed optimally
  template<typename PA>
  real_t ewaldEnergyLROpt_orig(const PA& R)
  {
    real_t ee = ewald_constant_lr_energy;

    real_t ve = 0.0;
    size_t n=0;
    for(const auto& k: kpoints)
    {
      real_t coskr = 0.0;
      real_t sinkr = 0.0;
      real_t c,s;
      size_t i=0;
      for(const auto& r: R)
      {
        real_t q = Q[i];
        sincos(dot(k,r),&s,&c);
        coskr += q*c;
        sinkr += q*s;
        i++;
      }
      real_t rhok2 = coskr*coskr + sinkr*sinkr;

      ve += rhok2*vk[n];
      n++;
    }

    ee += 0.5*ve - vk_sum;

    return ee;
  }




















  template<typename DT>
  real_t ewaldEnergySR0(const DT& dt)
  {
    const auto& dr = dt.getDisplacements();
    real_t ee = 0.0;
    RspaceEwaldPairPotential vsr(rexponent_ewald);
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
        ee += Q[i]*Q[j]*vsr(dr[i][j]);
    return ee;
  }






  template<typename DT>
  real_t ewaldEnergySRDT(const DT& dt)
  {
    auto& dr = dt.getDisplacements();
    real_t ee = 0.0;
    size_t Npairs = (N*(N-1))/2;
    IntVec nmax;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
      {
        real_t tol_ewald = tol/(Q[i]*Q[j])/Npairs;
        RspaceEwaldTerm rfunc(dr[i][j], A, rexponent_ewald);
        real_t rval = anisotropicGridSum(rfunc, nmax, true, tol_ewald);
        ee += Q[i]*Q[j]*rval;
      }
    return ee;
  }



  template<typename DT>
  real_t ewaldEnergyLRDT(const DT& dt)
  {
    auto& dr = dt.getDisplacements();
    real_t ee = ewald_constant_lr_energy;
    size_t Npairs = (N*(N-1))/2;
    IntVec nmax;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
      {
        real_t tol_ewald = tol/(Q[i]*Q[j])/Npairs;
        KspaceEwaldTerm kfunc(dr[i][j], B, kexponent_ewald, kprefactor_ewald);
        real_t kval = anisotropicGridSum(kfunc, nmax, false, tol_ewald);
        ee += Q[i]*Q[j]*kval;
      }
    return ee;
  }





};


} // namespace ewaldtools
} // namespace qmcplusplus
#endif
