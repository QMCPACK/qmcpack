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




/// General class to perform Ewald sums
class AnisotropicEwald
{
private:

  RealMat A;      // Real-space cell axes
  RealMat B;      // K-space cell axes
  ChargeArray Q;  // List of charges for all particles
  real_t tol;     // Tolerance in Ha for total energy adaptive sums

  size_t N;       // Number of particles
  real_t volume;  // Volume of real-space cell, \Omega in Drummond 2008
  real_t QQmax;   // Maximum charge product (used for setting tolerances)

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
  //std::rhok2                       // Charge structure factor

public:

  /// Empty constructor
  AnisotropicEwald() : nmax_anisotropic(0) { }

  /// State-initializing constructor
  AnisotropicEwald(const RealMat& A_in, const ChargeArray& Q_in, real_t tol_in = 1e-10,real_t kappa_in=-1.0)
    : A(A_in), Q(Q_in), tol(tol_in), nmax_anisotropic(0)
  {

    initialize(A_in,Q_in,tol_in,kappa_in);
  }

  /// Initialize constant data, including Madelung sums
  void initialize(const RealMat& A_in, const ChargeArray& Q_in, real_t tol_in = 1e-10,real_t kappa_in=-1.0)
  {
  
    A = A_in;
    Q = Q_in;
    tol = tol_in;

    // Reciprocal lattice vectors
    B = 2 * M_PI * transpose(inverse(A));

    // Set Ewald sum constants
    volume = std::abs(det(A));
    // Madelung constants
    if(kappa_in<0.0)
      kappa_madelung = getKappaMadelung(volume);
    else
      kappa_madelung = kappa_in;
    rexponent_madelung  = kappa_madelung;
    kexponent_madelung  = -1. / (4 * std::pow(kappa_madelung, 2));
    kprefactor_madelung = 4 * M_PI / volume;
    // Ewald constants
    if(kappa_in<0.0)
      kappa_ewald = getKappaEwald(volume);
    else
      kappa_ewald = kappa_in;
    rexponent_ewald  = 1. / (std::sqrt(2.) * kappa_ewald);
    kexponent_ewald  = -std::pow(kappa_ewald, 2) / 2;
    kprefactor_ewald = 4 * M_PI / volume;
    kconstant_ewald  = -2 * M_PI * std::pow(kappa_ewald, 2) / volume;

    // Number of particles
    N = Q.size();

    // Find maximum self-interaction charge product
    QQmax = 0.0;
    for (size_t i = 0; i < N; ++i)
      QQmax = std::max(std::abs(Q[i] * Q[i]), QQmax);

    // Calculate and store Madelung energy
    madelung_energy = madelungEnergy();

    // Calculate and store constant energy term for LR interaction
    real_t lrc = 0.0;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
        lrc += Q[i]*Q[j]*kconstant_ewald;
    ewald_constant_lr_energy = lrc;
  }


  /// Calculate the Madelung constant
  real_t madelungConstant() const
  {
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


  /// Update the current maximum anisotropic grid dimensions
  void updateNmax(const IntVec& nmax)
  {
    for(size_t d: {0,1,2})
      nmax_anisotropic[d] = std::max(nmax[d],nmax_anisotropic[d]);
  }


  /// Return the current maximum anisotropic grid dimensions
  IntVec getNmax()
  {
    return nmax_anisotropic;
  }


  /// Set the maximum anisotropic grid dimensions
  void setNmax(const IntVec& nmax)
  {
    nmax_anisotropic = nmax;
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
    updateNmax(nmax);
    // Compute the k-space sum (excludes zero)
    real_t kval = anisotropicGridSum(kfunc, nmax, false, tol_ewald);
    updateNmax(nmax);

    // Sum all contributions to get the Ewald pair interaction
    real_t epp = cval + rval + kval;

    return epp;
  }


  /// Total energy for all pairs using adaptive anisotropic grids
  template<typename PA>
  real_t ewaldEnergy(const PA& R)
  {
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
    real_t ee = 0.0;
    size_t Npairs = (N*(N-1))/2;
    IntVec nmax;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
      {
        real_t tol_ewald = tol/(Q[i]*Q[j])/Npairs;
        RspaceEwaldTerm rfunc(R[i]-R[j], A, rexponent_ewald);
        real_t rval = anisotropicGridSum(rfunc, nmax, true, tol_ewald);
        ee += Q[i]*Q[j]*rval;
      }
    return ee;
  }


  /// LR (k-space) part of total energy computed adaptively 
  template<typename PA>
  real_t ewaldEnergyLR(const PA& R)
  {
    real_t ee = ewald_constant_lr_energy;
    size_t Npairs = (N*(N-1))/2;
    IntVec nmax;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
      {
        real_t tol_ewald = tol/(Q[i]*Q[j])/Npairs;
        KspaceEwaldTerm kfunc(R[i]-R[j], B, kexponent_ewald, kprefactor_ewald);
        real_t kval = anisotropicGridSum(kfunc, nmax, false, tol_ewald);
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
  }


  void setupOpt()
  {
    setupOpt(nmax_anisotropic);
  }

  
 /// LR (k-space) part of total energy computed optimally
  template<typename PA>
  real_t ewaldEnergyLROpt(const PA& R)
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






















  template<typename PA>
  real_t ewaldEnergySR0(const PA& R)
  {
    real_t ee = 0.0;
    size_t Npairs = (N*(N-1))/2;
    IntVec nmax;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
      {
        real_t tol_ewald = tol/(Q[i]*Q[j])/Npairs;
        RspaceEwaldTerm rfunc(R[i]-R[j], A, rexponent_ewald);
        IntVec iv = 0;
        real_t rval = rfunc(iv);
        ee += Q[i]*Q[j]*rval;
      }
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



  template<typename DT>
  real_t ewaldEnergySR0DT(const DT& dt)
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
        IntVec iv = 0;
        real_t rval = rfunc(iv);
        ee += Q[i]*Q[j]*rval;
      }
    return ee;
  }





};


} // namespace ewaldtools
} // namespace qmcplusplus
#endif
