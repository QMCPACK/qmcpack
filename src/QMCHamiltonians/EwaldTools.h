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

  using namespace ewaldref;



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



class AnisotropicEwald
{
private:

  RealMat A;
  RealMat B;
  ChargeArray Q;
  real_t tol;

  size_t N;
  real_t volume;
  real_t QQmax;

  real_t kappa_madelung;
  real_t rconst_madelung;
  real_t kconst_madelung;
  real_t kfactor_madelung;

  real_t kappa_ewald;
  real_t rconst_ewald;
  real_t kconst_ewald;
  real_t kfactor_ewald;

  real_t madelung_constant;
  real_t madelung_energy;

  IntVec nmax_anisotropic;

public:

  AnisotropicEwald(const RealMat& A_in, const ChargeArray& Q_in, real_t tol_in = 1e-10,real_t kappa_in=-1.0)
    : A(A_in), Q(Q_in), tol(tol_in), nmax_anisotropic(0)
  {

    // Reciprocal lattice vectors
    B = 2 * M_PI * transpose(inverse(A));

    // Set Ewald sum constants
    volume = std::abs(det(A));
    // Madelung constants
    if(kappa_in<0.0)
      kappa_madelung = getKappaMadelung(volume);
    else
      kappa_madelung = kappa_in;
    rconst_madelung  = kappa_madelung;
    kconst_madelung  = -1. / (4 * std::pow(kappa_madelung, 2));
    kfactor_madelung = 4 * M_PI / volume;
    // Ewald constants
    if(kappa_in<0.0)
      kappa_ewald = getKappaEwald(volume);
    else
      kappa_ewald = kappa_in;
    rconst_ewald  = 1. / (std::sqrt(2.) * kappa_ewald);
    kconst_ewald  = -std::pow(kappa_ewald, 2) / 2;
    kfactor_ewald = 4 * M_PI / volume;

    // Number of particles
    N = Q.size();

    // Find maximum self-interaction charge product
    QQmax = 0.0;
    for (size_t i = 0; i < N; ++i)
      QQmax = std::max(std::abs(Q[i] * Q[i]), QQmax);

    // Calculate and store Madelung energy
    madelung_energy = madelungEnergy();
  }


  real_t madelungConstant()
  {
    // Set Madelung tolerance
    real_t tol_madelung = tol * 2. / QQmax;

    // Create real-/k-space fuctors for terms within the Madelung sums
    RspaceMadelungTerm rfunc(A, rconst_madelung);
    KspaceMadelungTerm kfunc(B, kconst_madelung, kfactor_madelung);
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


  void updateNmax(const IntVec& nmax)
  {
    for(size_t d: {0,1,2})
      nmax_anisotropic[d] = std::max(nmax[d],nmax_anisotropic[d]);
  }


  IntVec getNmax()
  {
    return nmax_anisotropic;
  }


  real_t ewaldPairPotential(const RealVec& r, real_t tol_ewald)
  {
    // Create real-/k-space fuctors for terms within the Ewald pair potential sums
    RspaceEwaldTerm rfunc(r, A, rconst_ewald);
    KspaceEwaldTerm kfunc(r, B, kconst_ewald, kfactor_ewald);
    IntVec nmax;

    // Compute the constant term
    real_t cval = -2 * M_PI * std::pow(kappa_ewald, 2) / volume;
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


  real_t ewaldEnergy(const PosArray& R)
  {
    real_t ee = madelung_energy;
    size_t Npairs = (N*(N-1))/2;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
        ee += Q[i]*Q[j]*ewaldPairPotential(R[i]-R[j],tol/(Q[i]*Q[j])/Npairs);
    return ee;
  }


  real_t fixedGridEwaldEnergy(const PosArray& R, const IntVec& nmax)
  {
    real_t ee = madelung_energy;
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < i; ++j)
      {
        RealVec r = R[i]-R[j];

        // Create real-/k-space fuctors for terms within the Ewald pair potential sums
        RspaceEwaldTerm rfunc(r, A, rconst_ewald);
        KspaceEwaldTerm kfunc(r, B, kconst_ewald, kfactor_ewald);
    
        // Compute the constant term
        real_t cval = -2 * M_PI * std::pow(kappa_ewald, 2) / volume;
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


  real_t fixedGridEwaldEnergy(const PosArray& R)
  {
    return fixedGridEwaldEnergy(R, nmax_anisotropic);
  }
};


} // namespace qmcplusplus
#endif
