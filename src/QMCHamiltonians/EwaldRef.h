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

#include <cmath>

#include "Configuration.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"
#include "Utilities/NewTimer.h"
#include "simd/allocator.hpp"


namespace qmcplusplus
{
namespace ewaldref
{
/// Reference Ewald implemented for 3D only
enum
{
  DIM = 3
};

/// Type for integers
using int_t = int;
/// Type for floating point numbers
using real_t = double;
/// Type for integer vectors of length DIM
using IntVec = TinyVector<int_t, DIM>;
/// Type for floating point vectors of length DIM
using RealVec = TinyVector<real_t, DIM>;
/// Type for floating point matrices of shape DIM,DIM
using RealMat = Tensor<real_t, DIM>;
/// Type for lists of particle positions
using PosArray = std::vector<RealVec>;
/// Type for lists of particle charges
using ChargeArray = std::vector<real_t>;


/// Functor for term within the real-space sum in Drummond 2008 formula 7
class RspaceMadelungTerm
{
private:
  /// The real-space cell axes
  const RealMat a;
  /// The constant \kappa in Drummond 2008 formula 7
  const real_t rconst;

public:
  RspaceMadelungTerm(const RealMat& a_in, real_t rconst_in) : a(a_in), rconst(rconst_in) {}

  real_t operator()(const IntVec& i) const
  {
    RealVec Rv = dot(i, a);
    real_t R   = std::sqrt(dot(Rv, Rv));
    real_t rm  = std::erfc(rconst * R) / R;
    return rm;
  }
};


/// Functor for term within the k-space sum in Drummond 2008 formula 7
class KspaceMadelungTerm
{
private:
  /// The k-space cell axes
  const RealMat b;
  /// The constant 1/(4\kappa^2) in Drummond 2008 formula 7
  const real_t kconst;
  /// The constant 4\pi/\Omega in Drummond 2008 formula 7
  const real_t kfactor;

public:
  KspaceMadelungTerm(const RealMat& b_in, real_t kconst_in, real_t kfactor_in)
      : b(b_in), kconst(kconst_in), kfactor(kfactor_in)
  {}

  real_t operator()(const IntVec& i) const
  {
    RealVec Kv = dot(i, b);
    real_t K2  = dot(Kv, Kv);
    real_t km  = kfactor * std::exp(kconst * K2) / K2;
    return km;
  }
};


/// Functor for term within the real-space sum in Drummond 2008 formula 6
class RspaceEwaldTerm
{
private:
  /// The inter-particle separation vector
  const RealVec r;
  /// The real-space cell axes
  const RealMat a;
  /// The constant 1/(\sqrt{2}kappa) in Drummond 2008 formula 6
  const real_t rconst;

public:
  RspaceEwaldTerm(RealVec r_in, RealMat a_in, real_t rconst_in) : r(r_in), a(a_in), rconst(rconst_in) {}

  real_t operator()(const IntVec& i) const
  {
    RealVec Rv = dot(i, a);
    for (int_t d : {0, 1, 2})
      Rv[d] -= r[d];
    real_t R  = std::sqrt(dot(Rv, Rv));
    real_t rm = std::erfc(rconst * R) / R;
    return rm;
  }
};

class RspaceEwaldTermForTile
{
private:
  const PosArray& R;
  const aligned_vector<real_t>& qqs;
  const size_t row_first;
  const size_t row_last;
  const size_t col_first;
  const size_t col_last;
  /// The inter-particle separation vector
  const RealVec r;
  /// The real-space cell axes
  const RealMat a;
  /// The constant 1/(\sqrt{2}kappa) in Drummond 2008 formula 6
  const real_t rconst;

public:
  RspaceEwaldTermForTile(const PosArray& R_in,
                  const aligned_vector<real_t>& qqs_in,
                  const size_t row_first_in,
                  const size_t row_last_in,
                  const size_t col_first_in,
                  const size_t col_last_in,
                  const RealMat& a_in,
                  const real_t rconst_in)
      : R(R_in),
        qqs(qqs_in),
        row_first(row_first_in),
        row_last(row_last_in),
        col_first(col_first_in),
        col_last(col_last_in),
        a(a_in),
        rconst(rconst_in) {}

  real_t getZeroContribution() const
  {
    real_t v(0);
    size_t icount(0);
    for (size_t i = row_first; i < row_last; ++i)
      #pragma omp simd
      for (size_t j = col_first; j < std::min(col_last, i); ++j)
      {
        RealVec r(R[i] - R[j]);
        real_t R  = std::sqrt(dot(r, r));
        real_t f = std::erfc(rconst * R) / R;
        v += f * qqs[icount];
        icount++;
      }
    return v;
  }

  real_t operator()(const IntVec& i, aligned_vector<real_t>& dva) const
  {
    real_t v(0);
    size_t icount(0);
    RealVec Rv = dot(i, a);
    for (size_t i = row_first; i < row_last; ++i)
      #pragma omp simd
      for (size_t j = col_first; j < std::min(col_last, i); ++j)
      {
        RealVec r(Rv - R[i] + R[j]);
        real_t R  = std::sqrt(dot(r, r));
        real_t f = std::erfc(rconst * R) / R;
        v += f * qqs[icount];
        dva[icount] += std::abs(f);
        icount++;
      }
    return v;
  }
};


/// Functor for term within the k-space sum in Drummond 2008 formula 6
class KspaceEwaldTerm
{
private:
  /// The inter-particle separation vector
  const RealVec r;
  /// The k-space cell axes
  const RealMat b;
  /// The constant -\kappa^2/2 in Drummond 2008 formula 6
  const real_t kconst;
  /// The constant 4\pi/\Omega in Drummond 2008 formula 6
  const real_t kfactor;

public:
  KspaceEwaldTerm(const RealVec& r_in, const RealMat& b_in, real_t kconst_in, real_t kfactor_in)
      : r(r_in), b(b_in), kconst(kconst_in), kfactor(kfactor_in)
  {}

  real_t operator()(const IntVec& i) const
  {
    RealVec Kv = dot(i, b);
    real_t K2  = dot(Kv, Kv);
    real_t Kr  = dot(Kv, r);
    real_t km  = kfactor * std::exp(kconst * K2) * std::cos(Kr) / K2;
    return km;
  }
};

inline size_t countPairs(const size_t row_first, const size_t row_last, const size_t col_first, const size_t col_last)
{
  size_t count(0);
  for (size_t i = row_first; i < row_last; ++i)
  {
    const size_t actual_col_last = std::min(col_last, i);
    if (col_first < actual_col_last)
      count += actual_col_last - col_first;
  }
  return count;
}

/// Functor implements r independent part of KspaceEwaldTerm
class KspaceEwaldTermForTile
{
private:
  const PosArray& R;
  const aligned_vector<real_t>& qqs;
  const size_t row_first;
  const size_t row_last;
  const size_t col_first;
  const size_t col_last;
  /// The k-space cell axes
  const RealMat b;
  /// The constant -\kappa^2/2 in Drummond 2008 formula 6
  const real_t kconst;
  /// The constant 4\pi/\Omega in Drummond 2008 formula 6
  const real_t kfactor;

  // scratch spaces
  aligned_vector<real_t> row_phase_sin;
  aligned_vector<real_t> row_phase_cos;
  aligned_vector<real_t> col_phase_sin;
  aligned_vector<real_t> col_phase_cos;
  aligned_vector<real_t> cosKvRij;

  real_t computeK2Exponetial(const IntVec& iv, RealVec& Kv) const
  {
    Kv        = dot(iv, b);
    real_t K2 = dot(Kv, Kv);
    return kfactor * std::exp(kconst * K2) / K2;
  }

public:
  /** constructor
   *  @param R Arary of particle positions
   *  @param qqs charge products of pairs
   *  @param row_first the first row particle index of a given tile
   *  @param row_last the last row particle index + 1 of a given tile
   *  @param col_first the first column particle index of a given tile
   *  @param col_last the last column particle index + 1 of a given tile
   *  @param b k-space cell axes.
   *  @param kconst The constant 1/(4\kappa^2) in Drummond 2008 formula 7
   *  @param kfactor The constant 4\pi/\Omega in Drummond 2008 formula 7
   */
  KspaceEwaldTermForTile(const PosArray& R_in,
                         const aligned_vector<real_t>& qqs_in,
                         const size_t row_first_in,
                         const size_t row_last_in,
                         const size_t col_first_in,
                         const size_t col_last_in,
                         const RealMat& b_in,
                         real_t kconst_in,
                         real_t kfactor_in)
      : R(R_in),
        qqs(qqs_in),
        row_first(row_first_in),
        row_last(row_last_in),
        col_first(col_first_in),
        col_last(col_last_in),
        b(b_in),
        kconst(kconst_in),
        kfactor(kfactor_in)
  {
    assert(row_last >= row_first);
    assert(col_last >= col_first);
    row_phase_sin.resize(row_last - row_first);
    row_phase_cos.resize(row_last - row_first);
    col_phase_sin.resize(col_last - col_first);
    col_phase_cos.resize(col_last - col_first);
    cosKvRij.resize(countPairs(row_first_in, row_last_in, col_first_in, col_last_in));
    assert(qqs.size() == cosKvRij.size());
  }

  real_t getZeroContribution() const { return 0; }

  real_t operator()(const IntVec& iv, aligned_vector<real_t>& dva)
  {
    assert(dva.size() == cosKvRij.size());
    auto* restrict qqs_ptr      = qqs.data();
    auto* restrict dva_ptr      = dva.data();
    auto* restrict cosKvRij_ptr = cosKvRij.data();

    auto* restrict row_cos = row_phase_cos.data();
    auto* restrict row_sin = row_phase_sin.data();
    auto* restrict col_cos = col_phase_cos.data();
    auto* restrict col_sin = col_phase_sin.data();

    RealVec Kv;
    real_t K2prefactor = computeK2Exponetial(iv, Kv);

#pragma omp simd aligned(row_cos, row_sin)
    for (size_t i = 0; i < row_last - row_first; ++i)
      sincos(dot(Kv, R[i + row_first]), &row_sin[i], &row_cos[i]);
#pragma omp simd aligned(col_cos, col_sin)
    for (size_t j = 0; j < col_last - col_first; ++j)
      sincos(-dot(Kv, R[j + col_first]), &col_sin[j], &col_cos[j]);

    real_t v(0);
    size_t icount = 0;
    for (size_t i = row_first; i < row_last; ++i)
#pragma omp simd aligned(row_cos, row_sin, col_cos, col_sin)
      for (size_t j = col_first; j < std::min(col_last, i); ++j)
      {
        cosKvRij[icount] =
            row_cos[i - row_first] * col_cos[j - col_first] - row_sin[i - row_first] * col_sin[j - col_first];
        icount++;
      }
#pragma omp simd aligned(cosKvRij_ptr, qqs_ptr, dva_ptr)
    for (size_t icount = 0; icount < cosKvRij.size(); ++icount)
    {
      real_t f = K2prefactor * cosKvRij_ptr[icount];
      v += f * qqs_ptr[icount];
      dva_ptr[icount] += std::abs(f);
    }
    return v;
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
real_t gridSum(T& function, bool zero = true, real_t tol = 1e-11)
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
  while (std::abs(dva) > tol)
  {
    dva = 0.0;
    dv  = 0.0; // Surface shell contribution.
    // Sum over new surface planes perpendicular to the x direction.
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
    // Sum over new surface planes perpendicular to the y direction.
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
    // Sum over new surface planes perpendicular to the z direction.
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
    v += dv;
  }

  return v;
}

/** vector check of individual tolerance
 */
template<typename T>
bool isLargerThanTol(const aligned_vector<T>& vals, const aligned_vector<T>& tols)
{
  for (size_t i = 0; i < vals.size(); i++)
    if (vals[i] > tols[i])
      return true;
  return false;
}

/** Similar to gridSum, perform a sum over successively larger cubic integer grid
 *  @param tol: Tolerance for the Ewald pair interaction in Ha.
 */
template<typename T>
real_t gridSumTile(T& function, const aligned_vector<real_t>& tols)
{
  aligned_vector<real_t> dva(tols.size(), std::numeric_limits<real_t>::max());

  real_t v = function.getZeroContribution();
  int_t indmax = 0;
  // Sum over cubic surface shells until the tolerance is reached.
  while (isLargerThanTol(dva, tols))
  {
    std::fill(dva.begin(), dva.end(), 0.0);

    // Surface shell contribution. The i,j,k loop can be threaded if the outer level threading is replaced with MPI.
    // Sum over new surface planes perpendicular to the x direction.
    for (int_t i : {-indmax - 1, indmax + 1})
      for (int_t j = -indmax; j < indmax + 1; ++j)
        for (int_t k = -indmax; k < indmax + 1; ++k)
          v += function(IntVec(i, j, k), dva);
    // Sum over new surface planes perpendicular to the y direction.
    for (int_t j : {-indmax - 1, indmax + 1})
      for (int_t k = -indmax; k < indmax + 1; ++k)
        for (int_t i = -indmax - 1; i < indmax + 2; ++i)
          v += function(IntVec(i, j, k), dva);
    // Sum over new surface planes perpendicular to the z direction.
    for (int_t k : {-indmax - 1, indmax + 1})
      for (int_t i = -indmax - 1; i < indmax + 2; ++i)
        for (int_t j = -indmax - 1; j < indmax + 2; ++j)
          v += function(IntVec(i, j, k), dva);
    indmax++;
  }

  return v;
}

real_t getKappaMadelung(real_t volume)
{
  real_t radius = std::pow(3. * volume / (4 * M_PI), 1. / 3);
  //return 2 * M_PI / ( 8 * radius );
  //return 1.0;
  return std::sqrt(M_PI) / radius;
}

/** Compute the Madelung constant to a given tolerance
 *
 *  Corresponds to the entirety of Drummond 2008 formula 7.
 *
 *  @param a: Real-space cell axes.
 *
 *  @param tol: Tolerance for the Madelung constant in Ha.
 */
real_t madelungSum(const RealMat& a, real_t tol = 1e-10)
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

  // Create real-/k-space fuctors for terms within the sums in formula 7
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
  //return 2 * M_PI / ( 8 * radius );
  //return 1.0;
  return radius / std::sqrt(2 * M_PI);
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
real_t ewaldSum(const RealVec& r, const RealMat& a, real_t tol = 1e-10)
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

  // Create real-/k-space fuctors for terms within the sums in formula 6
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

/** Similar to ewaldSum, compute the Ewald interaction of particle pairs in a tile
 *  @param a Real-space cell axes.
 *  @param R Arary of particle positions
 *  @param Q Arary of particle charges
 *  @param row_first the first row particle index of a given tile
 *  @param row_last the last row particle index + 1 of a given tile
 *  @param col_first the first column particle index of a given tile
 *  @param col_last the last column particle index + 1 of a given tile
 *  @param tol: Tolerance for the Ewald pair interaction in Ha.
 */
real_t ewaldSumTile(const RealMat& a,
                    const PosArray& R,
                    const ChargeArray& Q,
                    const size_t row_first,
                    const size_t row_last,
                    const size_t col_first,
                    const size_t col_last,
                    const real_t tol = 1e-10)
{
  // Real-space cell volume
  const real_t volume = std::abs(det(a));
  // k-space cell axes
  const RealMat b = 2 * M_PI * transpose(inverse(a));
  // k-space cutoff (kappa)
  real_t kconv = getKappaEwald(volume);

  // Set constants for real-/k-space Ewald pair functors
  const real_t rconst  = 1. / (std::sqrt(2.) * kconv);
  const real_t kconst  = -std::pow(kconv, 2) / 2;
  const real_t kfactor = 4 * M_PI / volume;

  size_t count = countPairs(row_first, row_last, col_first, col_last);

  aligned_vector<real_t> tols(count);
  aligned_vector<real_t> qqs(count);

  size_t icount = 0;
  for (size_t i = row_first; i < row_last; ++i)
    for (size_t j = col_first; j < std::min(col_last, i); ++j)
    {
      real_t qq    = Q[i] * Q[j];
      qqs[icount]  = qq;
      tols[icount] = tol / qq;
      icount++;
    }

  real_t es(0);
  { // Create k-space fuctors for terms within the sums in formula 6
    ScopedTimer KspaceTimer(TimerManager.createTimer("Kspace"));
    KspaceEwaldTermForTile kfuncTile(R, qqs, row_first, row_last, col_first, col_last, b, kconst, kfactor);
    es += gridSumTile(kfuncTile, tols);
  }
  { // Create real-space fuctors for terms within the sums in formula 7
    ScopedTimer RspaceTimer(TimerManager.createTimer("Rspace"));
    RspaceEwaldTermForTile rfuncTile(R, qqs, row_first, row_last, col_first, col_last, a, rconst);
    es += gridSumTile(rfuncTile, tols);
  }
  { // Compute the constant term
    const real_t cval = -2 * M_PI * kconv * kconv / volume;
    es += cval * accumulate(qqs.begin(), qqs.end(), real_t(0));
  }

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
real_t ewaldEnergy(const RealMat& a, const PosArray& R, const ChargeArray& Q, real_t tol = 1e-10)
{
  // Timer for EwaldRef
  ScopedTimer totalEwaldTimer(TimerManager.createTimer("EwaldRef"));

  // Number of particles
  const size_t N = R.size();

  // Total Ewald potential energy
  real_t ve = 0.0;

  {
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

  if (false)
  {
    ScopedTimer EwaldSumTimer(TimerManager.createTimer("EwaldSum"));
    real_t ve_ewald(0);
// Sum the interaction terms for all particle pairs
#pragma omp parallel for reduction(+ : ve_ewald)
    for (size_t i = 1; i < N / 2 + 1; ++i)
    {
      for (size_t j = 0; j < i; ++j)
      {
        const real_t qq = Q[i] * Q[j];
        ve_ewald += qq * ewaldSum(R[i] - R[j], a, tol / qq);
      }

      const size_t i_reverse = N - i;
      if (i == i_reverse)
        continue;

      for (size_t j = 0; j < i_reverse; ++j)
      {
        const real_t qq = Q[i_reverse] * Q[j];
        ve_ewald += qq * ewaldSum(R[i_reverse] - R[j], a, tol / qq);
      }
    }
    ve += ve_ewald;
  }
  else
  {
    ScopedTimer EwaldSumTimer(TimerManager.createTimer("EwaldSumOpt"));

    size_t tile_size, max_n_tile_1d, n_tiles;
    if (N <= 1)
    {
      // protect N-2 underflow
      tile_size     = 1;
      max_n_tile_1d = 0;
      n_tiles       = 0;
    }
    else
    {
      const size_t n_threads = omp_get_max_threads();
      // making a guess to obtain intial tile_size as large as possible.
      tile_size = (N - 2) / static_cast<size_t>(std::ceil(std::sqrt(n_threads * 2 + 0.25) - 0.5));
      // searching for the largest possible n_tiles <= n_threads
      do
      {
        tile_size++;
        max_n_tile_1d = (N - 2) / tile_size + 1;
        n_tiles       = max_n_tile_1d * (max_n_tile_1d + 1) / 2;
      } while (n_tiles > n_threads);
    }

    /* The computation of the lower triangle of N^2 matrix ( actually N(N-1) pairs) are sliced into square tiles of tile_size^2.
     * Larger tiles are preferred for maximizing data reuse but a limited number of tiles are bad for parallelization.
     * The above code searches for the largest possible number of tiles to fulfill parallel execution unit
     * and also make tile_size as large as possible.
     * Currently tiles are distributed over threads. If needed, we can replace threads with MPI ranks at this level.
     * Distributing tiles in the lower triangle may not be the optimal strategy for load balancing. It can be improved if needed.
     */
    std::vector<TinyVector<size_t, 4>> work_list;
    work_list.reserve(n_tiles);
    for (size_t i = 1; i < N; i += tile_size)
      for (size_t j = 0; j < i; j += tile_size)
        work_list.emplace_back(i, std::min(i + tile_size, N), j, std::min(j + tile_size, N - 1));

    real_t ve_ewald(0);
#pragma omp parallel for reduction(+ : ve_ewald)
    for (size_t i = 0; i < work_list.size(); i++)
      ve_ewald += ewaldSumTile(a, R, Q, work_list[i][0], work_list[i][1], work_list[i][2], work_list[i][3], tol);
    ve += ve_ewald;
  }

  return ve;
}

} // namespace ewaldref
} // namespace qmcplusplus

#endif
