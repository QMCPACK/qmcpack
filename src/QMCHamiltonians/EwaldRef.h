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


namespace qmcplusplus
{

/// Rewrap the OHMMS_DIM constant as DIM
enum
{
  DIM = OHMMS_DIM
};

/// Type for integers
using IntType     = int;
/// Type for floating point numbers
using RealType    = QMCTraits::RealType;
/// Type for integer vectors of length DIM
using IntVec      = TinyVector<IntType, DIM>;
/// Type for floating point vectors of length DIM
using RealVec     = TinyVector<RealType, DIM>;
/// Type for floating point matrices of shape DIM,DIM
using RealMat     = Tensor<RealType, DIM>;
/// Type for lists of particle positions
using PosArray    = std::vector<RealVec>;
/// Type for lists of particle charges
using ChargeArray = std::vector<RealType>;


/// Functor for term within the real-space sum in Drummond 2008 formula 7
class RspaceMadelungTerm
{
  public:

  /// The real-space cell axes
  RealMat a;
  /// The constant \kappa in Drummond 2008 formula 7
  RealType rconst;

  RspaceMadelungTerm(RealMat a_in,RealType rconst_in)
  {
    a      = a_in;
    rconst = rconst_in;
  }

  RealType operator()(IntVec i)
  {
    RealVec Rv  = dot(i,a);
    RealType R  = std::sqrt(dot(Rv,Rv));
    RealType rm = std::erfc(rconst*R)/R;
    return rm;
  }
};


/// Functor for term within the k-space sum in Drummond 2008 formula 7
class KspaceMadelungTerm
{
  public:

  /// The k-space cell axes
  RealMat b;
  /// The constant 1/(4\kappa^2) in Drummond 2008 formula 7
  RealType kconst;
  /// The constant 4\pi/\Omega in Drummond 2008 formula 7
  RealType kfactor;

  KspaceMadelungTerm(RealMat b_in,RealType kconst_in,RealType kfactor_in)
  {
    b       = b_in;
    kconst  = kconst_in;
    kfactor = kfactor_in;
  }

  RealType operator()(IntVec i)
  {
    RealVec  Kv = dot(i,b);
    RealType K2 = dot(Kv,Kv);
    RealType km = kfactor*std::exp(kconst*K2)/K2;
    return km;
  }
};


/// Functor for term within the real-space sum in Drummond 2008 formula 6
class RspaceEwaldTerm
{
  public:

  /// The inter-particle separation vector
  RealVec r;
  /// The real-space cell axes
  RealMat a;
  /// The constant 1/(\sqrt{2}kappa) in Drummond 2008 formula 6
  RealType rconst;

  RspaceEwaldTerm(RealVec r_in,RealMat a_in,RealType rconst_in)
  {
    r      = r_in;
    a      = a_in;
    rconst = rconst_in;
  }

  RealType operator()(IntVec i)
  {
    RealVec Rv = dot(i,a);
    for(IntType d: {0,1,2})
      Rv[d] -= r[d];
    RealType R  = std::sqrt(dot(Rv,Rv));
    RealType rm = std::erfc(rconst*R)/R;
    return rm;
  }
};


/// Functor for term within the k-space sum in Drummond 2008 formula 6
class KspaceEwaldTerm
{
  public:

  /// The inter-particle separation vector
  RealVec r;
  /// The k-space cell axes
  RealMat b;
  /// The constant -\kappa^2/2 in Drummond 2008 formula 6
  RealType kconst;
  /// The constant 4\pi/\Omega in Drummond 2008 formula 6
  RealType kfactor;

  KspaceEwaldTerm(RealVec r_in,RealMat b_in,RealType kconst_in,RealType kfactor_in)
  {
    r       = r_in;
    b       = b_in;
    kconst  = kconst_in;
    kfactor = kfactor_in;
  }

  RealType operator()(IntVec i)
  {
     RealVec  Kv = dot(i,b);
     RealType K2 = dot(Kv,Kv);
     RealType Kr = dot(Kv,r);
     RealType km = kfactor*std::exp(kconst*K2)*std::cos(Kr)/K2;
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
RealType gridSum(T& function,bool zero=true,RealType tol=1e-11)
{
  RealType dv = 1e99;
  RealType v  = 0.0;
  IntType im  = 0;
  IntType jm  = 0;
  IntType km  = 0;
  IntVec  iv;
  // Add the value at the origin, if requested.
  if(zero)
  {
    iv = 0;
    v += function(iv);
  }
  // Sum over cubic surface shells until the tolerance is reached.
  while(std::abs(dv)>tol)
  {
    dv = 0.0; // Surface shell contribution.
    // Sum over new surface planes perpendicular to the x direction.
    im += 1;
    for(IntType i: {-im,im})
      for(IntType j=-jm; j<jm+1; ++j)
        for(IntType k=-km; k<km+1; ++k)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          dv += function(iv);
        }
    // Sum over new surface planes perpendicular to the y direction.
    jm += 1;
    for(IntType j: {-jm,jm})
      for(IntType k=-km; k<km+1; ++k)
        for(IntType i=-im; i<im+1; ++i)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          dv += function(iv);
        }
    // Sum over new surface planes perpendicular to the z direction.
    km += 1;
    for(IntType k: {-km,km})
      for(IntType i=-im; i<im+1; ++i)
        for(IntType j=-jm; j<jm+1; ++j)
        {
          iv[0] = i;
          iv[1] = j;
          iv[2] = k;
          dv += function(iv);
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
RealType madelungSum(RealMat a,RealType tol=1e-10)
{
  // Real-space cell volume
  RealType volume = std::abs(det(a));
  // k-space cell axes
  RealMat b       = 2*M_PI*transpose(inverse(a));
  // real space cutoff
  RealType rconv  = 8*std::pow(3.*volume/(4*M_PI),1./3);
  // k-space cutoff (kappa)
  RealType kconv  = 2*M_PI/rconv;

  // Set constants for real-/k-space Madelung functors
  RealType rconst  = kconv;
  RealType kconst  = -1./(4*std::pow(kconv,2));
  RealType kfactor = 4*M_PI/volume;

  // Create real-/k-space fuctors for terms within the sums in formula 7
  RspaceMadelungTerm rfunc(a,rconst);
  KspaceMadelungTerm kfunc(b,kconst,kfactor);

  // Compute the constant term
  RealType cval = -M_PI/(std::pow(kconv,2)*volume)-2*kconv/std::sqrt(M_PI);
  // Compute the real-space sum (excludes zero)
  RealType rval = gridSum(rfunc,false,tol);
  // Compute the k-space sum (excludes zero)
  RealType kval = gridSum(kfunc,false,tol);

  // Sum all contributions to get the Madelung constant
  RealType ms = cval + rval + kval;

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
RealType ewaldSum(RealVec r,RealMat a,RealType tol=1e-10)
{
  // Real-space cell volume
  RealType volume = std::abs(det(a));
  // k-space cell axes
  RealMat b       = 2*M_PI*transpose(inverse(a));
  // real space cutoff
  RealType rconv  = 8*std::pow(3.*volume/(4*M_PI),1./3);
  // k-space cutoff (kappa)
  RealType kconv  = 2*M_PI/rconv;

  // Set constants for real-/k-space Ewald pair functors
  RealType rconst  = 1./(std::sqrt(2.)*kconv);
  RealType kconst  = -std::pow(kconv,2)/2;
  RealType kfactor = 4*M_PI/volume;

  // Create real-/k-space fuctors for terms within the sums in formula 6
  RspaceEwaldTerm rfunc(r,a,rconst);
  KspaceEwaldTerm kfunc(r,b,kconst,kfactor);

  // Compute the constant term
  RealType cval = -2*M_PI*std::pow(kconv,2)/volume;
  // Compute the real-space sum (includes zero)
  RealType rval = gridSum(rfunc,true,tol);
  // Compute the k-space sum (excludes zero)
  RealType kval = gridSum(kfunc,false,tol);

  // Sum all contributions to get the Ewald pair interaction
  RealType es = cval + rval + kval;

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
RealType ewaldEnergy(RealMat a,PosArray R,ChargeArray Q,RealType tol=1e-10)
{
  // Number of particles
  IntType N = R.size();

  // Maximum self-interaction charge product
  RealType qqmax=0.0;
  for(size_t i=0; i<N; ++i)
    qqmax = std::max(std::abs(Q[i]*Q[i]),qqmax);

  // Compute the Madelung term (Drummond 2008 formula 7)
  RealType vm = madelungSum(a,tol*2./qqmax);

  // Total Ewald potential energy
  RealType ve = 0.0;

  // Sum the Madelung self interaction for each particle
  for(size_t i=0; i<N; ++i)
    ve += Q[i]*Q[i]*vm/2;

  // Sum the interaction terms for all particle pairs
  for(size_t i=0; i<N; ++i)
    for(size_t j=0; j<i; ++j)
    {
      RealType qq = Q[i]*Q[j];
      ve += qq*ewaldSum(R[i]-R[j],a,tol/qq);
    }

  return ve;
}

}

#endif
