//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_BSPLINE_FUNCTOR_H
#define QMCPLUSPLUS_BSPLINE_FUNCTOR_H

#include <cstdio>
#include <memory>
#include "OptimizableFunctorBase.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Numerics/LinearFit.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

namespace qmcplusplus
{

/**BsplineFunctor class for the Jastrows
 * REAL is the real type used by offload target, it is the correct type for the mw data pointers
 * and is also used to coerce/implicitly convert the Real type inherited OptimizableFunctorBase into that buffer
 * if offload is off this happens too but is just an implementation quirk.
 */
template<typename REAL>
struct BsplineFunctor : public OptimizableFunctorBase
{
  using Real = OptimizableFunctorBase::real_type;

  static constexpr Real A0 = -1.0 / 6.0, A1 = 3.0 / 6.0, A2 = -3.0 / 6.0, A3 = 1.0 / 6.0;
  static constexpr Real A4 = 3.0 / 6.0, A5 = -6.0 / 6.0, A6 = 0.0 / 6.0, A7 = 4.0 / 6.0;
  static constexpr Real A8 = -3.0 / 6.0, A9 = 3.0 / 6.0, A10 = 3.0 / 6.0, A11 = 1.0 / 6.0;
  static constexpr Real A12 = 1.0 / 6.0, A13 = 0.0 / 6.0, A14 = 0.0 / 6.0, A15 = 0.0 / 6.0;

  static constexpr Real dA0 = 0.0, dA1 = -0.5, dA2 = 1.0, dA3 = -0.5;
  static constexpr Real dA4 = 0.0, dA5 = 1.5, dA6 = -2.0, dA7 = 0.0;
  static constexpr Real dA8 = 0.0, dA9 = -1.5, dA10 = 1.0, dA11 = 0.5;
  static constexpr Real dA12 = 0.0, dA13 = 0.5, dA14 = 0.0, dA15 = 0.0;

  static constexpr Real d2A0 = 0.0, d2A1 = 0.0, d2A2 = -1.0, d2A3 = 1.0;
  static constexpr Real d2A4 = 0.0, d2A5 = 0.0, d2A6 = 3.0, d2A7 = -2.0;
  static constexpr Real d2A8 = 0.0, d2A9 = 0.0, d2A10 = -3.0, d2A11 = 1.0;
  static constexpr Real d2A12 = 0.0, d2A13 = 0.0, d2A14 = 1.0, d2A15 = 0.0;

  static constexpr Real d3A0 = 0.0, d3A1 = 0.0, d3A2 = 0.0, d3A3 = -1.0;
  static constexpr Real d3A4 = 0.0, d3A5 = 0.0, d3A6 = 0.0, d3A7 = 3.0;
  static constexpr Real d3A8 = 0.0, d3A9 = 0.0, d3A10 = 0.0, d3A11 = -3.0;
  static constexpr Real d3A12 = 0.0, d3A13 = 0.0, d3A14 = 0.0, d3A15 = 1.0;

  std::shared_ptr<Vector<Real, OffloadAllocator<Real>>> spline_coefs_;

  int NumParams;
  Real DeltaR, DeltaRInv;
  Real CuspValue;
  Real Y, dY, d2Y;
  // Stores the derivatives w.r.t. coefs
  // of the u, du/dr, and d2u/dr2
  std::vector<TinyVector<Real, 3>> SplineDerivs;
  std::vector<Real> Parameters;
  std::vector<std::string> ParameterNames;
  std::string elementType, pairType;
  std::string fileName;

  bool notOpt;
  bool periodic;

  ///constructor
  BsplineFunctor(const std::string& my_name, Real cusp = 0.0)
      : OptimizableFunctorBase(my_name), NumParams(0), CuspValue(cusp), notOpt(false), periodic(true)
  {
    cutoff_radius = 0.0;
  }

  OptimizableFunctorBase* makeClone() const override { return new BsplineFunctor(*this); }

  void setCusp(Real c) override { CuspValue = c; }

  void setPeriodic(bool p) override { periodic = p; }

  void resize(int n)
  {
    NumParams    = n;
    int numCoefs = NumParams + 4;
    int numKnots = numCoefs - 2;
    DeltaR       = cutoff_radius / (Real)(numKnots - 1);
    DeltaRInv    = 1.0 / DeltaR;
    Parameters.resize(n);
    spline_coefs_ = std::make_shared<Vector<Real, OffloadAllocator<Real>>>(numCoefs);
    SplineDerivs.resize(numCoefs);
  }

  /** reset coefs from Parameters
   */
  void reset() override
  {
    const int numCoefs = NumParams + 4;
    const int numKnots = numCoefs - 2;
    DeltaR             = cutoff_radius / (Real)(numKnots - 1);
    DeltaRInv          = 1.0 / DeltaR;
    auto& coefs        = *spline_coefs_;
    for (int i = 0; i < coefs.size(); i++)
      coefs[i] = 0.0;
    // Ensure that cusp conditions is satisfied at the origin
    coefs[1] = Parameters[0];
    coefs[2] = Parameters[1];
    coefs[0] = Parameters[1] - 2.0 * DeltaR * CuspValue;
    for (int i = 2; i < Parameters.size(); i++)
      coefs[i + 1] = Parameters[i];
    coefs.updateTo();
  }

  /** compute value, first and second derivatives for [iStart, iEnd) pairs
   * @param iat the source particle that should be avoided (self pairs)
   * @param iStart starting particle index
   * @param iEnd ending particle index
   * @param _distArray distance arrUay
   * @param _valArray  u(r_j) for j=[iStart,iEnd)
   * @param _gradArray  du(r_j)/dr /r_j for j=[iStart,iEnd)
   * @param _lapArray  d2u(r_j)/dr2 for j=[iStart,iEnd)
   * @param distArrayCompressed temp storage to filter r_j < cutoff_radius
   * @param distIndices temp storage for the compressed index
   */
  void evaluateVGL(const int iat,
                   const int iStart,
                   const int iEnd,
                   const REAL* _distArray,
                   REAL* restrict _valArray,
                   REAL* restrict _gradArray,
                   REAL* restrict _laplArray,
                   REAL* restrict distArrayCompressed,
                   int* restrict distIndices) const;

  /** compute value, gradient and laplacian for target particles
   * This more than just a batched call of evaluateVGL
   * @param iat the source particle that should be avoided (self pairs)
   * @param num_groups the number of source particle groups
   * @param functors for the num_groups of source particles
   * @param n_src the number of source particles
   * @param grp_ids the group ids of the n_src source particles
   * @param nw batch size (number of walkers)
   * @param mw_vgl return resutls. Multi walker value, gradient and laplacian [nw][1(v)+DIM(g)+1(l)]
   * @param n_padded the padded size of source particles
   * @param mw_dist Multi walker distance table [nw][1(distance)+DIM(displacements)][n_padded]
   * @param mw_cur_allu Multi walker value, first and second derivatives of pair potentials [nw][DIM][n_padded]. if mw_cur_allu is dual space, only update device side.
   * @param transfer_buffer temporary transfer buffer.
   *
   * If mw_dist is dual space, up-to-date data is assumed on device.
   * If mw_cur_allu is dual space, data is created on the device and there is no transfer to the host
   * because it will be consumed by mw_updateVGL on the device.
   */
  static void mw_evaluateVGL(const int iat,
                             const int num_groups,
                             const BsplineFunctor* const functors[],
                             const int n_src,
                             const int* grp_ids,
                             const int nw,
                             REAL* mw_vgl, // [nw][DIM+2]
                             const int n_padded,
                             const REAL* mw_dist, // [nw][DIM+1][n_padded]
                             REAL* mw_cur_allu,   // [nw][3][n_padded]
                             Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer);

  /** evaluate sum of the pair potentials for [iStart,iEnd)
   * @param iat dummy
   * @param iStart starting particle index
   * @param iEnd ending particle index
   * @param _distArray distance arrUay
   * @param distArrayCompressed temp storage to filter r_j < cutoff_radius
   * @return \f$\sum u(r_j)\f$ for r_j < cutoff_radius
   */
  REAL evaluateV(const int iat,
                 const int iStart,
                 const int iEnd,
                 const REAL* restrict _distArray,
                 REAL* restrict distArrayCompressed) const;

  /** compute value for target-source particle pair potentials
   * This more than just a batched call of evaluateV
   * @param num_groups the number of source particle groups
   * @param functors for the num_groups of source particles
   * @param n_src the number of source particles
   * @param grp_ids the group ids of the n_src source particles
   * @param nnum_pairs the number of particle pairs
   * @param ref_at the source particles that should be avoided (self pairs)
   * @param mw_vgl return resutls. Multi walker value, gradient and laplacian [nw][1(v)+DIM(g)+1(l)]
   * @param dist_stride the offset of distance pointers between to consecutive walkers
   * @param mw_dist Multi walker distance table [nw][1(distance)+DIM(displacements)][n_padded]
   * @param transfer_buffer temporary transfer buffer.
   *
   * If mw_dist is dual space, up-to-date data is assumed on device.
   */
  static void mw_evaluateV(const int num_groups,
                           const BsplineFunctor* const functors[],
                           const int n_src,
                           const int* grp_ids,
                           const int num_pairs,
                           const int* ref_at,
                           const REAL* mw_dist,
                           const int dist_stride,
                           REAL* mw_vals,
                           Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer);

  inline static Real evaluate_impl(Real r, const Real* coefs, const Real DeltaRInv)
  {
    r *= DeltaRInv;
    REAL ipart;
    const REAL t = std::modf(r, &ipart);
    const int i  = (int)ipart;

    Real sCoef0 = coefs[i + 0];
    Real sCoef1 = coefs[i + 1];
    Real sCoef2 = coefs[i + 2];
    Real sCoef3 = coefs[i + 3];

    return (sCoef0 * (((A0 * t + A1) * t + A2) * t + A3) + sCoef1 * (((A4 * t + A5) * t + A6) * t + A7) +
            sCoef2 * (((A8 * t + A9) * t + A10) * t + A11) + sCoef3 * (((A12 * t + A13) * t + A14) * t + A15));
  }

  inline Real evaluate(Real r) const
  {
    Real u(0);
    if (r < cutoff_radius)
      u = evaluate_impl(r, spline_coefs_->data(), DeltaRInv);
    return u;
  }

  inline Real evaluate(Real r, Real rinv) { return Y = evaluate(r, dY, d2Y); }

  inline void evaluateAll(Real r, Real rinv) { Y = evaluate(r, dY, d2Y); }

  inline static Real evaluate_impl(Real r, const Real* coefs, const Real DeltaRInv, Real& dudr, Real& d2udr2)
  {
    r *= DeltaRInv;
    REAL ipart;
    const REAL t = std::modf(r, &ipart);
    const int i  = (int)ipart;

    Real sCoef0 = coefs[i + 0];
    Real sCoef1 = coefs[i + 1];
    Real sCoef2 = coefs[i + 2];
    Real sCoef3 = coefs[i + 3];

    d2udr2 = DeltaRInv * DeltaRInv *
        (sCoef0 * (d2A2 * t + d2A3) + sCoef1 * (d2A6 * t + d2A7) + sCoef2 * (d2A10 * t + d2A11) +
         sCoef3 * (d2A14 * t + d2A15));

    dudr = DeltaRInv *
        (sCoef0 * ((dA1 * t + dA2) * t + dA3) + sCoef1 * ((dA5 * t + dA6) * t + dA7) +
         sCoef2 * ((dA9 * t + dA10) * t + dA11) + sCoef3 * ((dA13 * t + dA14) * t + dA15));

    Real u = (sCoef0 * (((A0 * t + A1) * t + A2) * t + A3) + sCoef1 * (((A4 * t + A5) * t + A6) * t + A7) +
              sCoef2 * (((A8 * t + A9) * t + A10) * t + A11) + sCoef3 * (((A12 * t + A13) * t + A14) * t + A15));
    return u;
  }

  inline Real evaluate(Real r, Real& dudr, Real& d2udr2)
  {
    Real u(0);
    dudr   = Real(0);
    d2udr2 = Real(0);

    if (r < cutoff_radius)
      u = evaluate_impl(r, spline_coefs_->data(), DeltaRInv, dudr, d2udr2);
    return u;
  }


  inline Real evaluate(Real r, Real& dudr, Real& d2udr2, Real& d3udr3)
  {
    if (r >= cutoff_radius)
    {
      dudr = d2udr2 = d3udr3 = 0.0;
      return 0.0;
    }
    // Real eps = 1.0e-5;
    //       Real dudr_FD = (evaluate(r+eps)-evaluate(r-eps))/(2.0*eps);
    //       Real d2udr2_FD = (evaluate(r+eps)+evaluate(r-eps)-2.0*evaluate(r))/(eps*eps);
    // Real d3udr3_FD = (-1.0*evaluate(r+1.0*eps)
    //         +2.0*evaluate(r+0.5*eps)
    //         -2.0*evaluate(r-0.5*eps)
    //         +1.0*evaluate(r-1.0*eps))/(eps*eps*eps);
    r *= DeltaRInv;
    Real ipart, t;
    t     = std::modf(r, &ipart);
    int i = (int)ipart;
    Real tp[4];
    tp[0]       = t * t * t;
    tp[1]       = t * t;
    tp[2]       = t;
    tp[3]       = 1.0;
    auto& coefs = *spline_coefs_;
    d3udr3      = DeltaRInv * DeltaRInv * DeltaRInv *
        (coefs[i + 0] * (d3A0 * tp[0] + d3A1 * tp[1] + d3A2 * tp[2] + d3A3 * tp[3]) +
         coefs[i + 1] * (d3A4 * tp[0] + d3A5 * tp[1] + d3A6 * tp[2] + d3A7 * tp[3]) +
         coefs[i + 2] * (d3A8 * tp[0] + d3A9 * tp[1] + d3A10 * tp[2] + d3A11 * tp[3]) +
         coefs[i + 3] * (d3A12 * tp[0] + d3A13 * tp[1] + d3A14 * tp[2] + d3A15 * tp[3]));
    d2udr2 = DeltaRInv * DeltaRInv *
        (coefs[i + 0] * (d2A0 * tp[0] + d2A1 * tp[1] + d2A2 * tp[2] + d2A3 * tp[3]) +
         coefs[i + 1] * (d2A4 * tp[0] + d2A5 * tp[1] + d2A6 * tp[2] + d2A7 * tp[3]) +
         coefs[i + 2] * (d2A8 * tp[0] + d2A9 * tp[1] + d2A10 * tp[2] + d2A11 * tp[3]) +
         coefs[i + 3] * (d2A12 * tp[0] + d2A13 * tp[1] + d2A14 * tp[2] + d2A15 * tp[3]));
    dudr = DeltaRInv *
        (coefs[i + 0] * (dA0 * tp[0] + dA1 * tp[1] + dA2 * tp[2] + dA3 * tp[3]) +
         coefs[i + 1] * (dA4 * tp[0] + dA5 * tp[1] + dA6 * tp[2] + dA7 * tp[3]) +
         coefs[i + 2] * (dA8 * tp[0] + dA9 * tp[1] + dA10 * tp[2] + dA11 * tp[3]) +
         coefs[i + 3] * (dA12 * tp[0] + dA13 * tp[1] + dA14 * tp[2] + dA15 * tp[3]));
    //       if (std::abs(dudr_FD-dudr) > 1.0e-8)
    //  std::cerr << "Error in BsplineFunction:  dudr = " << dudr
    //       << "  dudr_FD = " << dudr_FD << std::endl;
    //       if (std::abs(d2udr2_FD-d2udr2) > 1.0e-4)
    //  std::cerr << "Error in BsplineFunction:  r = " << r << "  d2udr2 = " << dudr
    //       << "  d2udr2_FD = " << d2udr2_FD << "  rcut = " << cutoff_radius << std::endl;
    // if (std::abs(d3udr3_FD-d3udr3) > 1.0e-4)
    //  std::cerr << "Error in BsplineFunction:  r = " << r << "  d3udr3 = " << dudr
    //       << "  d3udr3_FD = " << d3udr3_FD << "  rcut = " << cutoff_radius << std::endl;
    return (coefs[i + 0] * (A0 * tp[0] + A1 * tp[1] + A2 * tp[2] + A3 * tp[3]) +
            coefs[i + 1] * (A4 * tp[0] + A5 * tp[1] + A6 * tp[2] + A7 * tp[3]) +
            coefs[i + 2] * (A8 * tp[0] + A9 * tp[1] + A10 * tp[2] + A11 * tp[3]) +
            coefs[i + 3] * (A12 * tp[0] + A13 * tp[1] + A14 * tp[2] + A15 * tp[3]));
  }

  /** update value, gradient and laplacian for target particles
   * It serves multile walkers and handles update in a batched fashion
   * @param iat the source particle that should be avoided (self pairs)
   * @param isAccepted accept/reject status
   * @param num_groups the number of source particle groups
   * @param functors for the num_groups of source particles
   * @param n_src the number of source particles
   * @param grp_ids the group ids of the n_src source particles
   * @param nw batch size (number of walkers)
   * @param mw_vgl Multi walker value, gradient and laplacian [nw][1(v)+DIM(g)+1(l)]
   * @param n_padded the padded size of source particles
   * @param mw_dist Multi walker distance table [new + old][nw][1(distance)+DIM(displacements)][n_padded]
   * @param mw_allUat, returned results. Multi walker value, gradient and laplacian of pair potentials [nw][1(v)+DIM(g)+1(l)][n_padded]
   * @param mw_cur_allu Multi walker value, first and second derivatives of pair potentials [nw][DIM][n_padded]
   * @param transfer_buffer temporary transfer buffer
   *
   * If mw_dist is dual space, up-to-date data is assumed on device.
   * If mw_cur_allu is dual space, data on the device is consumed and no transfer is needed.
   */
  static void mw_updateVGL(const int iat,
                           const std::vector<bool>& isAccepted,
                           const int num_groups,
                           const BsplineFunctor* const functors[],
                           const int n_src,
                           const int* grp_ids,
                           const int nw,
                           REAL* mw_vgl, // [nw][DIM+2]
                           const int n_padded,
                           const REAL* mw_dist, // [nw][DIM+1][n_padded]
                           REAL* mw_allUat,     // [nw][DIM+2][n_padded]
                           REAL* mw_cur_allu,   // [nw][3][n_padded]
                           Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer);

  inline bool evaluateDerivatives(Real r, std::vector<TinyVector<Real, 3>>& derivs) override
  {
    if (r >= cutoff_radius)
      return false;
    r *= DeltaRInv;
    Real ipart, t;
    t     = std::modf(r, &ipart);
    int i = (int)ipart;
    Real tp[4];
    tp[0] = t * t * t;
    tp[1] = t * t;
    tp[2] = t;
    tp[3] = 1.0;

    auto& coefs     = *spline_coefs_;
    SplineDerivs[0] = TinyVector<Real, 3>(0.0);
    // d/dp_i u(r)
    SplineDerivs[i + 0][0] = A0 * tp[0] + A1 * tp[1] + A2 * tp[2] + A3 * tp[3];
    SplineDerivs[i + 1][0] = A4 * tp[0] + A5 * tp[1] + A6 * tp[2] + A7 * tp[3];
    SplineDerivs[i + 2][0] = A8 * tp[0] + A9 * tp[1] + A10 * tp[2] + A11 * tp[3];
    SplineDerivs[i + 3][0] = A12 * tp[0] + A13 * tp[1] + A14 * tp[2] + A15 * tp[3];
    // d/dp_i du/dr
    SplineDerivs[i + 0][1] = DeltaRInv * (dA1 * tp[1] + dA2 * tp[2] + dA3 * tp[3]);
    SplineDerivs[i + 1][1] = DeltaRInv * (dA5 * tp[1] + dA6 * tp[2] + dA7 * tp[3]);
    SplineDerivs[i + 2][1] = DeltaRInv * (dA9 * tp[1] + dA10 * tp[2] + dA11 * tp[3]);
    SplineDerivs[i + 3][1] = DeltaRInv * (dA13 * tp[1] + dA14 * tp[2] + dA15 * tp[3]);
    // d/dp_i d2u/dr2
    SplineDerivs[i + 0][2] = DeltaRInv * DeltaRInv * (d2A2 * tp[2] + d2A3 * tp[3]);
    SplineDerivs[i + 1][2] = DeltaRInv * DeltaRInv * (d2A6 * tp[2] + d2A7 * tp[3]);
    SplineDerivs[i + 2][2] = DeltaRInv * DeltaRInv * (d2A10 * tp[2] + d2A11 * tp[3]);
    SplineDerivs[i + 3][2] = DeltaRInv * DeltaRInv * (d2A14 * tp[2] + d2A15 * tp[3]);

    int imin = std::max(i, 1);
    int imax = std::min(i + 4, NumParams + 1);
    for (int n = imin; n < imax; ++n)
      derivs[n - 1] = SplineDerivs[n];
    derivs[1] += SplineDerivs[0];

    //Real v[4],dv[4],d2v[4];
    //v[0] = A[ 0]*tp[0] + A[ 1]*tp[1] + A[ 2]*tp[2] + A[ 3]*tp[3];
    //v[1] = A[ 4]*tp[0] + A[ 5]*tp[1] + A[ 6]*tp[2] + A[ 7]*tp[3];
    //v[2] = A[ 8]*tp[0] + A[ 9]*tp[1] + A10*tp[2] + A11*tp[3];
    //v[3] = A12*tp[0] + A13*tp[1] + A14*tp[2] + A15*tp[3];
    //// d/dp_i du/dr
    //dv[0] = DeltaRInv * (dA[ 1]*tp[1] + dA[ 2]*tp[2] + dA[ 3]*tp[3]);
    //dv[1] = DeltaRInv * (dA[ 5]*tp[1] + dA[ 6]*tp[2] + dA[ 7]*tp[3]);
    //dv[2] = DeltaRInv * (dA[ 9]*tp[1] + dA10*tp[2] + dA11*tp[3]);
    //dv[3] = DeltaRInv * (dA13*tp[1] + dA14*tp[2] + dA15*tp[3]);
    //// d/dp_i d2u/dr2
    //d2v[0] = DeltaRInv * DeltaRInv * (d2A[ 2]*tp[2] + d2A[ 3]*tp[3]);
    //d2v[1] = DeltaRInv * DeltaRInv * (d2A[ 6]*tp[2] + d2A[ 7]*tp[3]);
    //d2v[2] = DeltaRInv * DeltaRInv * (d2A10*tp[2] + d2A11*tp[3]);
    //d2v[3] = DeltaRInv * DeltaRInv * (d2A14*tp[2] + d2A15*tp[3]);

    //int imin=std::max(i,1);
    //int imax=std::min(i+4,NumParams+1)-1;
    //int n=imin-1, j=imin-i;
    //while(n<imax && j<4)
    //{
    //  derivs[n] = TinyVector<Real,3>(v[j],dv[j],d2v[j]);
    //  n++; j++;
    //}
    //if(i==0) derivs[1]+= TinyVector<Real,3>(v[0],dv[0],d2v[0]);

    return true;
  }

  inline bool evaluateDerivatives(Real r, std::vector<Real>& derivs) override
  {
    if (r >= cutoff_radius)
      return false;
    Real tp[4], v[4], ipart, t;
    t        = std::modf(r * DeltaRInv, &ipart);
    tp[0]    = t * t * t;
    tp[1]    = t * t;
    tp[2]    = t;
    tp[3]    = 1.0;
    v[0]     = A0 * tp[0] + A1 * tp[1] + A2 * tp[2] + A3 * tp[3];
    v[1]     = A4 * tp[0] + A5 * tp[1] + A6 * tp[2] + A7 * tp[3];
    v[2]     = A8 * tp[0] + A9 * tp[1] + A10 * tp[2] + A11 * tp[3];
    v[3]     = A12 * tp[0] + A13 * tp[1] + A14 * tp[2] + A15 * tp[3];
    int i    = (int)ipart;
    int imin = std::max(i, 1);
    int imax = std::min(i + 4, NumParams + 1) - 1;
    int n = imin - 1, j = imin - i;
    while (n < imax && j < 4)
    {
      derivs[n] = v[j];
      n++;
      j++;
    }
    if (i == 0)
      derivs[1] += v[0];
    return true;
  }

  inline Real f(Real r) override
  {
    if (r >= cutoff_radius)
      return 0.0;
    return evaluate(r);
  }
  inline Real df(Real r) override
  {
    if (r >= cutoff_radius)
      return 0.0;
    Real du, d2u;
    evaluate(r, du, d2u);
    return du;
  }


  bool put(xmlNodePtr cur) override
  {
    ReportEngine PRE("BsplineFunctor", "put(xmlNodePtr)");
    //CuspValue = -1.0e10;
    NumParams = 0;
    //cutoff_radius = 0.0;
    OhmmsAttributeSet rAttrib;
    Real radius = -1.0;
    rAttrib.add(NumParams, "size");
    rAttrib.add(radius, "rcut");
    rAttrib.add(radius, "cutoff");
    rAttrib.put(cur);
    if (radius < 0.0)
      if (periodic)
      {
        app_log() << "    Jastrow cutoff unspecified.  Setting to Wigner-Seitz radius = " << cutoff_radius << std::endl;
        app_log() << std::endl;
      }
      else
      {
        APP_ABORT("  Jastrow cutoff unspecified.  Cutoff must be given when using open boundary conditions");
      }
    else if (periodic && radius > cutoff_radius)
    {
      if (radius - cutoff_radius > 1e-4)
      {
        APP_ABORT("  The Jastrow cutoff specified should not be larger than Wigner-Seitz radius.");
      }
      else
      {
        app_log() << "  The Jastrow cutoff specified is slightly larger than the Wigner-Seitz radius.";
        app_log() << "  Setting to Wigner-Seitz radius = " << cutoff_radius << ".\n";
      }
    }
    else
      cutoff_radius = radius;
    if (NumParams == 0)
    {
      PRE.error("You must specify a positive number of parameters for the Bspline jastrow function.", true);
    }
    app_summary() << "     Number of parameters: " << NumParams << std::endl;
    app_summary() << "     Cusp: " << CuspValue << std::endl;
    app_summary() << "     Cutoff radius: " << cutoff_radius << std::endl;
    resize(NumParams);
    // Now read coefficents
    xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
    while (xmlCoefs != NULL)
    {
      std::string cname((const char*)xmlCoefs->name);
      if (cname == "coefficients")
      {
        std::string type("0"), id("0");
        std::string optimize("yes");
        OhmmsAttributeSet cAttrib;
        cAttrib.add(id, "id");
        cAttrib.add(type, "type");
        cAttrib.add(optimize, "optimize");
        cAttrib.put(xmlCoefs);
        if (type != "Array")
        {
          PRE.error("Unknown correlation type " + type + " in BsplineFunctor." + "Resetting to \"Array\"");
          xmlNewProp(xmlCoefs, (const xmlChar*)"type", (const xmlChar*)"Array");
        }
        std::vector<Real> params;
        putContent(params, xmlCoefs);
        if (params.size() == NumParams)
          Parameters = params;
        else
        {
          app_log() << "    Changing number of Bspline parameters from " << params.size() << " to " << NumParams
                    << ".  Performing fit:\n";
          // Fit function to new number of parameters
          const int numPoints = 500;
          BsplineFunctor<REAL> tmp_func("tmp_func", CuspValue);
          tmp_func.cutoff_radius = cutoff_radius;
          tmp_func.resize(params.size());
          tmp_func.Parameters = params;
          tmp_func.reset();
          std::vector<Real> y(numPoints);
          Matrix<Real> basis(numPoints, NumParams);
          std::vector<TinyVector<Real, 3>> derivs(NumParams);
          for (int i = 0; i < numPoints; i++)
          {
            Real r = (Real)i / (Real)numPoints * cutoff_radius;
            y[i]   = tmp_func.evaluate(r);
            evaluateDerivatives(r, derivs);
            for (int j = 0; j < NumParams; j++)
              basis(i, j) = derivs[j][0];
          }
          resize(NumParams);
          LinearFit(y, basis, Parameters);
          app_log() << "New parameters are:\n";
          for (int i = 0; i < Parameters.size(); i++)
            app_log() << "   " << Parameters[i] << std::endl;
        }
        if (optimize == "yes")
        {
          notOpt = false;
        }
        else
        {
          notOpt = true;
        }
        for (int i = 0; i < NumParams; i++)
        {
          std::stringstream sstr;
          sstr << id << "_" << i;
          myVars.insert(sstr.str(), (Real)Parameters[i], !notOpt, optimize::LOGLINEAR_P);
        }
        int left_pad_space = 5;
        app_log() << std::endl;
        myVars.print(app_log(), left_pad_space, true);
      }
      xmlCoefs = xmlCoefs->next;
    }
    reset();
    Real zeros = 0;
    for (int i = 0; i < NumParams; i++)
      zeros += Parameters[i] * Parameters[i];
    return zeros > 1.0e-12; //true if Parameters are not zero
  }

  void initialize(int numPoints,
                  std::vector<Real>& x,
                  std::vector<Real>& y,
                  REAL cusp,
                  REAL rcut,
                  std::string& id,
                  std::string& optimize)
  {
    ReportEngine PRE("BsplineFunctor", "initialize");
    NumParams     = numPoints;
    cutoff_radius = rcut;
    CuspValue     = cusp;
    if (NumParams == 0)
    {
      PRE.error("You must specify a positive number of parameters for the Bspline jastrow function.", true);
    }
    app_log() << "Initializing BsplineFunctor from array. \n";
    app_log() << " size = " << NumParams << " parameters " << std::endl;
    app_log() << " cusp = " << CuspValue << std::endl;
    app_log() << " rcut = " << cutoff_radius << std::endl;
    resize(NumParams);
    int npts = x.size();
    Matrix<Real> basis(npts, NumParams);
    std::vector<TinyVector<Real, 3>> derivs(NumParams);
    for (int i = 0; i < npts; i++)
    {
      Real r = x[i];
      if (r > cutoff_radius)
      {
        PRE.error("Error in BsplineFunctor::initialize: r > cutoff_radius.", true);
      }
      evaluateDerivatives(r, derivs);
      for (int j = 0; j < NumParams; j++)
        basis(i, j) = derivs[j][0];
    }
    resize(NumParams);
    LinearFit(y, basis, Parameters);
    app_log() << "New parameters are:\n";
    for (int i = 0; i < Parameters.size(); i++)
      app_log() << "   " << Parameters[i] << std::endl;
#if !defined(QMC_BUILD_SANDBOX_ONLY)
    if (optimize == "yes")
    {
      // Setup parameter names
      for (int i = 0; i < NumParams; i++)
      {
        std::stringstream sstr;
        sstr << id << "_" << i;
        myVars.insert(sstr.str(), (Real)Parameters[i], true, optimize::LOGLINEAR_P);
      }
      myVars.print(app_log());
    }
    else
#endif
    {
      notOpt = true;
      app_log() << "Parameters of BsplineFunctor id:" << id << " are not being optimized.\n";
    }
    reset();
  }

  void reportStatus(std::ostream& os) override
  {
    if (notOpt)
      return;
    myVars.print(os);
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    if (notOpt)
      return;
    myVars.getIndex(active);
  }

  void checkInVariablesExclusive(opt_variables_type& active) override
  {
    if (notOpt)
      return;
    myVars.setIndexDefault();
    active.insertFrom(myVars);
  }

  void resetParametersExclusive(const opt_variables_type& active) override
  {
    if (notOpt)
      return;
    for (int i = 0; i < Parameters.size(); ++i)
    {
      int loc = myVars.where(i);
      if (loc >= 0)
        Parameters[i] = std::real(myVars[i] = active[loc]);
    }
    reset();
  }

  // check if this object has active optimizable parameters
  bool isOptimizable()
  {
    if (notOpt)
      return false;
    for (int i = 0; i < Parameters.size(); ++i)
    {
      int loc = myVars.where(i);
      if (loc >= 0)
        return true;
    }
    return false;
  }
};

template<typename REAL>
inline REAL BsplineFunctor<REAL>::evaluateV(const int iat,
                                            const int iStart,
                                            const int iEnd,
                                            const REAL* restrict _distArray,
                                            REAL* restrict distArrayCompressed) const
{
  const Real* restrict distArray = _distArray + iStart;

  ASSUME_ALIGNED(distArrayCompressed);
  int iCount       = 0;
  const int iLimit = iEnd - iStart;

#pragma vector always
  for (int jat = 0; jat < iLimit; jat++)
  {
    Real r = distArray[jat];
    // pick the distances smaller than the cutoff and avoid the reference atom
    if (r < cutoff_radius && iStart + jat != iat)
      distArrayCompressed[iCount++] = distArray[jat];
  }

  Real d      = 0.0;
  auto& coefs = *spline_coefs_;
#pragma omp simd reduction(+ : d)
  for (int jat = 0; jat < iCount; jat++)
  {
    Real r = distArrayCompressed[jat];
    r *= DeltaRInv;
    const int i  = (int)r;
    const Real t = r - Real(i);
    Real d1      = coefs[i + 0] * (((A0 * t + A1) * t + A2) * t + A3);
    Real d2      = coefs[i + 1] * (((A4 * t + A5) * t + A6) * t + A7);
    Real d3      = coefs[i + 2] * (((A8 * t + A9) * t + A10) * t + A11);
    Real d4      = coefs[i + 3] * (((A12 * t + A13) * t + A14) * t + A15);
    d += (d1 + d2 + d3 + d4);
  }
  return d;
}

template<typename REAL>
inline void BsplineFunctor<REAL>::evaluateVGL(const int iat,
                                              const int iStart,
                                              const int iEnd,
                                              const REAL* _distArray,
                                              REAL* restrict _valArray,
                                              REAL* restrict _gradArray,
                                              REAL* restrict _laplArray,
                                              REAL* restrict distArrayCompressed,
                                              int* restrict distIndices) const
{
  Real dSquareDeltaRinv = DeltaRInv * DeltaRInv;
  constexpr Real cOne(1);

  //    START_MARK_FIRST();

  ASSUME_ALIGNED(distIndices);
  ASSUME_ALIGNED(distArrayCompressed);
  int iCount            = 0;
  int iLimit            = iEnd - iStart;
  const REAL* distArray = _distArray + iStart;
  REAL* valArray        = _valArray + iStart;
  REAL* gradArray       = _gradArray + iStart;
  REAL* laplArray       = _laplArray + iStart;

#pragma vector always
  for (int jat = 0; jat < iLimit; jat++)
  {
    Real r = distArray[jat];
    if (r < cutoff_radius && iStart + jat != iat)
    {
      distIndices[iCount]         = jat;
      distArrayCompressed[iCount] = r;
      iCount++;
    }
  }

  auto& coefs = *spline_coefs_;
#pragma omp simd
  for (int j = 0; j < iCount; j++)
  {
    Real r       = distArrayCompressed[j];
    int iScatter = distIndices[j];
    Real rinv    = cOne / r;
    r *= DeltaRInv;
    const int iGather = (int)r;
    const Real t      = r - Real(iGather);

    Real sCoef0 = coefs[iGather + 0];
    Real sCoef1 = coefs[iGather + 1];
    Real sCoef2 = coefs[iGather + 2];
    Real sCoef3 = coefs[iGather + 3];

    laplArray[iScatter] = dSquareDeltaRinv *
        (sCoef0 * (d2A2 * t + d2A3) + sCoef1 * (d2A6 * t + d2A7) + sCoef2 * (d2A10 * t + d2A11) +
         sCoef3 * (d2A14 * t + d2A15));

    gradArray[iScatter] = DeltaRInv * rinv *
        (sCoef0 * ((dA1 * t + dA2) * t + dA3) + sCoef1 * ((dA5 * t + dA6) * t + dA7) +
         sCoef2 * ((dA9 * t + dA10) * t + dA11) + sCoef3 * ((dA13 * t + dA14) * t + dA15));

    valArray[iScatter] =
        (sCoef0 * (((A0 * t + A1) * t + A2) * t + A3) + sCoef1 * (((A4 * t + A5) * t + A6) * t + A7) +
         sCoef2 * (((A8 * t + A9) * t + A10) * t + A11) + sCoef3 * (((A12 * t + A13) * t + A14) * t + A15));
  }
}

extern template struct BsplineFunctor<QMCTraits::RealType>;

} // namespace qmcplusplus
#endif
