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


#include "QMCHamiltonians/EwaldRef.h"


namespace qmcplusplus
{
namespace ewaldref
{


real_t getKappaMadelung(real_t volume)
{
  real_t radius = std::pow(3. * volume / (4 * M_PI), 1. / 3);
  return std::sqrt(M_PI) / radius;
}


real_t madelungSum(const RealMat& a, real_t tol)
{
  // Real-space cell volume
  real_t volume = std::abs(det(a));
  // k-space cell axes
  RealMat b = 2 * M_PI * transpose(inverse(a));
  // k-space cutoff (kappa)
  real_t kconv = getKappaMadelung(volume);

  // Set constants for real-/k-space Madelung functors
  real_t rconst  = kconv;
  real_t kconst  = -1. / (4 * std::pow(kconv, 2));
  real_t kfactor = 4 * M_PI / volume;

  // Create real-/k-space functors for terms within the sums in formula 7
  RspaceMadelungTerm rfunc(a, rconst);
  KspaceMadelungTerm kfunc(b, kconst, kfactor);

  // Compute the constant term
  real_t cval = -M_PI / (std::pow(kconv, 2) * volume) - 2 * kconv / std::sqrt(M_PI);
  // Compute the real-space sum (excludes zero)
  real_t rval = gridSum(rfunc, false, tol);
  // Compute the k-space sum (excludes zero)
  real_t kval = gridSum(kfunc, false, tol);

  // Sum all contributions to get the Madelung constant
  real_t ms = cval + rval + kval;

  return ms;
}


real_t getKappaEwald(real_t volume)
{
  real_t radius = std::pow(3. * volume / (4 * M_PI), 1. / 3);
  return radius / std::sqrt(2 * M_PI);
}


real_t ewaldSum(const RealVec& r, const RealMat& a, real_t tol)
{
  // Real-space cell volume
  real_t volume = std::abs(det(a));
  // k-space cell axes
  RealMat b = 2 * M_PI * transpose(inverse(a));
  // k-space cutoff (kappa)
  real_t kconv = getKappaEwald(volume);

  // Set constants for real-/k-space Ewald pair functors
  real_t rconst  = 1. / (std::sqrt(2.) * kconv);
  real_t kconst  = -std::pow(kconv, 2) / 2;
  real_t kfactor = 4 * M_PI / volume;

  // Create real-/k-space functors for terms within the sums in formula 6
  RspaceEwaldTerm rfunc(r, a, rconst);
  KspaceEwaldTerm kfunc(r, b, kconst, kfactor);

  // Compute the constant term
  real_t cval = -2 * M_PI * std::pow(kconv, 2) / volume;
  // Compute the real-space sum (includes zero)
  real_t rval = gridSum(rfunc, true, tol);
  // Compute the k-space sum (excludes zero)
  real_t kval = gridSum(kfunc, false, tol);

  // Sum all contributions to get the Ewald pair interaction
  real_t es = cval + rval + kval;

  return es;
}


real_t ewaldEnergy(const RealMat& a, const PosArray& R, const ChargeArray& Q, real_t tol)
{
  // Timer for EwaldRef
  ScopedTimer totalEwaldTimer(TimerManager.createTimer("EwaldRef"));

  // Number of particles
  const size_t N = R.size();

  // Total Ewald potential energy
  real_t ve = 0.0;

  {
    // Sum Madelung contributions
    ScopedTimer totalMadelungTimer(TimerManager.createTimer("MadelungSum"));
    // Maximum self-interaction charge product
    real_t qqmax = 0.0;
    for (size_t i = 0; i < N; ++i)
      qqmax = std::max(std::abs(Q[i] * Q[i]), qqmax);

    // Compute the Madelung term (Drummond 2008 formula 7)
    real_t vm = madelungSum(a, tol * 2. / qqmax);

    // Sum the Madelung self interaction for each particle
    for (size_t i = 0; i < N; ++i)
      ve += Q[i] * Q[i] * vm / 2;
  }

  {
    // Sum the interaction terms for all particle pairs
    ScopedTimer EwaldSumTimer(TimerManager.createTimer("EwaldSum"));

    int_t Npairs = (N*(N-1))/2;

    std::vector<real_t>  qq(Npairs);
    for(size_t i=0,n=0; i<N; ++i)
      for(size_t j=0; j<i; ++j,++n)
        qq[n] = Q[i]*Q[j];

    std::vector<RealVec> rr(Npairs);
    for(size_t i=0,n=0; i<N; ++i)
      for(size_t j=0; j<i; ++j,++n)
        rr[n] = R[i]-R[j];

#pragma omp parallel for reduction(+ : ve)
    for(size_t n=0; n<Npairs; ++n)
      ve += qq[n]*ewaldSum(rr[n],a,tol/qq[n]/Npairs);
  }

  return ve;
}


}
}
