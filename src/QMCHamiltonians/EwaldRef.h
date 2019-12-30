//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

/**@file EwaldRef.h
 *
 * @brief Computes Ewald sums of the potential energy to a given 
 *    tolerance for arbitrary collections of charges.
 *
 * The implementation follows formulas 6 and 7 from:
 *
 *   N. D. Drummond et al., Physical Review B 78 125106 (2008)
 *
 *   DOI:  https://doi.org/10.1103/PhysRevB.78.125106
 */

#ifndef EWALD_REF_H
#define EWALD_REF_H

#include<cmath>

#include "Configuration.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "Utilities/NewTimer.h"


namespace qmcplusplus
{

/// Reference Ewald implemented for 3D only
enum
{
  DIM = 3
};

/// Type for integers
using int_t     = int;
/// Type for floating point numbers
using real_t    = double;
/// Type for integer vectors of length DIM
using IntVec      = TinyVector<int_t, DIM>;
/// Type for floating point vectors of length DIM
using RealVec     = TinyVector<real_t, DIM>;
/// Type for floating point matrices of shape DIM,DIM
using RealMat     = Tensor<real_t, DIM>;
/// Type for lists of particle positions
using PosArray    = std::vector<RealVec>;
/// Type for lists of particle charges
using ChargeArray = std::vector<real_t>;


/// Functor for term within the real-space sum in Drummond 2008 formula 7
class RspaceMadelungTerm
{
  private:

  /// The real-space cell axes
  RealMat a;
  /// The constant \kappa in Drummond 2008 formula 7
  real_t rconst;

  public:

  RspaceMadelungTerm(RealMat a_in,real_t rconst_in)
  {
    a      = a_in;
    rconst = rconst_in;
  }

  real_t operator()(IntVec i)
  {
    RealVec Rv  = dot(i,a);
    real_t R  = std::sqrt(dot(Rv,Rv));
    real_t rm = std::erfc(rconst*R)/R;
    return rm;
  }
};


/// Functor for term within the k-space sum in Drummond 2008 formula 7
class KspaceMadelungTerm
{
  private:

  /// The k-space cell axes
  RealMat b;
  /// The constant 1/(4\kappa^2) in Drummond 2008 formula 7
  real_t kconst;
  /// The constant 4\pi/\Omega in Drummond 2008 formula 7
  real_t kfactor;

  public:

  KspaceMadelungTerm(RealMat b_in,real_t kconst_in,real_t kfactor_in)
  {
    b       = b_in;
    kconst  = kconst_in;
    kfactor = kfactor_in;
  }

  real_t operator()(IntVec i)
  {
    RealVec  Kv = dot(i,b);
    real_t K2 = dot(Kv,Kv);
    real_t km = kfactor*std::exp(kconst*K2)/K2;
    return km;
  }
};


/// Functor for term within the real-space sum in Drummond 2008 formula 6
class RspaceEwaldTerm
{
  private:

  /// The inter-particle separation vector
  RealVec r;
  /// The real-space cell axes
  RealMat a;
  /// The constant 1/(\sqrt{2}kappa) in Drummond 2008 formula 6
  real_t rconst;

  public:

  RspaceEwaldTerm(RealVec r_in,RealMat a_in,real_t rconst_in)
  {
    r      = r_in;
    a      = a_in;
    rconst = rconst_in;
  }

  real_t operator()(IntVec i)
  {
    RealVec Rv = dot(i,a);
    for(int_t d: {0,1,2})
      Rv[d] -= r[d];
    real_t R  = std::sqrt(dot(Rv,Rv));
    real_t rm = std::erfc(rconst*R)/R;
    return rm;
  }
};


/// Functor for term within the k-space sum in Drummond 2008 formula 6
class KspaceEwaldTerm
{
  private:

  /// The inter-particle separation vector
  RealVec r;
  /// The k-space cell axes
  RealMat b;
  /// The constant -\kappa^2/2 in Drummond 2008 formula 6
  real_t kconst;
  /// The constant 4\pi/\Omega in Drummond 2008 formula 6
  real_t kfactor;

  public:

  KspaceEwaldTerm(RealVec r_in,RealMat b_in,real_t kconst_in,real_t kfactor_in)
  {
    r       = r_in;
    b       = b_in;
    kconst  = kconst_in;
    kfactor = kfactor_in;
  }

  real_t operator()(IntVec i)
  {
     RealVec  Kv = dot(i,b);
     real_t K2 = dot(Kv,Kv);
     real_t Kr = dot(Kv,r);
     real_t km = kfactor*std::exp(kconst*K2)*std::cos(Kr)/K2;
     return km;
  }
};


/** Perform a sum over successively larger cubic integer grids
 *  in DIM dimensional space for arbitrary functors.
 *
 *  @param function: A functor accepting a point in the grid and 
 *    returning the real-valued contribution to the sum from that 
 *    point.
 *
 *  @param zero: Include the origin in the sum (or not).
 *
 *  @param tol: Tolerance for the sum.  Summation ceases when the 
 *    contribution to the sum from the outermost cubic shell of 
 *    points is less than tol.
 */
template<typename T>
real_t gridSum(T& function,bool zero=true,real_t tol=1e-11)
{
  real_t dv  = std::numeric_limits<real_t>::max();
  real_t dva = std::numeric_limits<real_t>::max();
  real_t v   = 0.0;
  int_t im   = 0;
  int_t jm   = 0;
  int_t km   = 0;
  IntVec  iv;
  // Add the value at the origin, if requested.
  if(zero)
  {
    iv = 0;
    v += function(iv);
  }
  // Sum over cubic surface shells until the tolerance is reached.
  while(std::abs(dva)>tol)
  {
    dva = 0.0;
    dv  = 0.0; // Surface shell contribution.
    // Sum over new surface planes perpendicular to the x direction.
    im += 1;
    for(int_t i: {-im,im})
      for(int_t j=-jm; j<jm+1; ++j)
        for(int_t k=-km; k<km+1; ++k)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          real_t f = function(iv);
          dv  += f;
          dva += std::abs(f);
        }
    // Sum over new surface planes perpendicular to the y direction.
    jm += 1;
    for(int_t j: {-jm,jm})
      for(int_t k=-km; k<km+1; ++k)
        for(int_t i=-im; i<im+1; ++i)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          real_t f = function(iv);
          dv  += f;
          dva += std::abs(f);
        }
    // Sum over new surface planes perpendicular to the z direction.
    km += 1;
    for(int_t k: {-km,km})
      for(int_t i=-im; i<im+1; ++i)
        for(int_t j=-jm; j<jm+1; ++j)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          real_t f = function(iv);
          dv  += f;
          dva += std::abs(f);
        }
    v += dv;
  }

  return v;
}


/** Compute the Madelung constant to a given tolerance
 *
 *  Corresponds to the entirety of Drummond 2008 formula 7.
 *
 *  @param a: Real-space cell axes.
 *
 *  @param tol: Tolerance for the Madelung constant in Ha.
 */
real_t madelungSum(RealMat a,real_t tol=1e-10)
{
  // Real-space cell volume
  real_t volume = std::abs(det(a));
  // k-space cell axes
  RealMat b       = 2*M_PI*transpose(inverse(a));
  // real space cutoff
  real_t rconv  = 2*M_PI;
  // k-space cutoff (kappa)
  real_t kconv  = 2*M_PI/rconv;

  // Set constants for real-/k-space Madelung functors
  real_t rconst  = kconv;
  real_t kconst  = -1./(4*std::pow(kconv,2));
  real_t kfactor = 4*M_PI/volume;

  // Create real-/k-space fuctors for terms within the sums in formula 7
  RspaceMadelungTerm rfunc(a,rconst);
  KspaceMadelungTerm kfunc(b,kconst,kfactor);

  // Compute the constant term
  real_t cval = -M_PI/(std::pow(kconv,2)*volume)-2*kconv/std::sqrt(M_PI);
  // Compute the real-space sum (excludes zero)
  real_t rval = gridSum(rfunc,false,tol);
  // Compute the k-space sum (excludes zero)
  real_t kval = gridSum(kfunc,false,tol);

  // Sum all contributions to get the Madelung constant
  real_t ms = cval + rval + kval;

  return ms;
}


/** Compute the Ewald interaction of a particle pair to a given tolerance
 *
 *  Corresponds to the entirety of Drummond 2008 formula 6.
 *
 *  @param r: Inter-particle separation vector.
 *
 *  @param a: Real-space cell axes.
 *
 *  @param tol: Tolerance for the Ewald pair interaction in Ha.
 */
real_t ewaldSum(RealVec r,RealMat a,real_t tol=1e-10)
{
  // Real-space cell volume
  real_t volume = std::abs(det(a));
  // k-space cell axes
  RealMat b       = 2*M_PI*transpose(inverse(a));
  // real space cutoff
  real_t rconv  = 2*M_PI;
  // k-space cutoff (kappa)
  real_t kconv  = 2*M_PI/rconv;

  // Set constants for real-/k-space Ewald pair functors
  real_t rconst  = 1./(std::sqrt(2.)*kconv);
  real_t kconst  = -std::pow(kconv,2)/2;
  real_t kfactor = 4*M_PI/volume;

  // Create real-/k-space fuctors for terms within the sums in formula 6
  RspaceEwaldTerm rfunc(r,a,rconst);
  KspaceEwaldTerm kfunc(r,b,kconst,kfactor);

  // Compute the constant term
  real_t cval = -2*M_PI*std::pow(kconv,2)/volume;
  // Compute the real-space sum (includes zero)
  real_t rval = gridSum(rfunc,true,tol);
  // Compute the k-space sum (excludes zero)
  real_t kval = gridSum(kfunc,false,tol);

  // Sum all contributions to get the Ewald pair interaction
  real_t es = cval + rval + kval;

  return es;
}


/** Compute the total Ewald potential energy for a collection of charges
 *
 *  Corresponds to the entirety of Drummond 2008 formula 5, but for 
 *    arbitrary charges.
 *
 *  @param a: Real-space cell axes.
 *
 *  @param R: List of particle coordinates.
 *
 *  @param R: List of particle charges.
 *
 *  @param tol: Tolerance for the total potential energy in Ha.
 */
real_t ewaldEnergy(RealMat a,PosArray R,ChargeArray Q,real_t tol=1e-10)
{
  // Timer for EwaldRef
  NewTimer& EwaldTimer(*TimerManager.createTimer("EwaldRef"));
  ScopedTimer totalEwaldTimer(&EwaldTimer);

  // Number of particles
  int_t N = R.size();

  // Maximum self-interaction charge product
  real_t qqmax=0.0;
  for(size_t i=0; i<N; ++i)
    qqmax = std::max(std::abs(Q[i]*Q[i]),qqmax);

  // Compute the Madelung term (Drummond 2008 formula 7)
  real_t vm = madelungSum(a,tol*2./qqmax);

  // Total Ewald potential energy
  real_t ve = 0.0;

  // Sum the Madelung self interaction for each particle
  for(size_t i=0; i<N; ++i)
    ve += Q[i]*Q[i]*vm/2;

  // Sum the interaction terms for all particle pairs
  for(size_t i=0; i<N; ++i)
    for(size_t j=0; j<i; ++j)
    {
      real_t qq = Q[i]*Q[j];
      ve += qq*ewaldSum(R[i]-R[j],a,tol/qq);
    }

  return ve;
}

} // namespace qmcplusplus

#endif
