//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *  Assume that coeffs.D1 and the LogLightGrid
 *   r_values.size() are equal
 *  Therefore r must be < r_max
 */

#ifndef QMCPLUSPLUS_MULTI_FUNCTOR_QUINTIC_SPLINE_SET_H
#define QMCPLUSPLUS_MULTI_FUNCTOR_QUINTIC_SPLINE_SET_H

#include <algorithm>

#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimQuinticSpline.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include <cassert>

namespace qmcplusplus
{
template<typename T>
struct LogGridLight
{
  T lower_bound;
  T upper_bound;
  T OneOverLogDelta;
  double LogDelta;
  std::vector<T> r_values;

  inline void set(T ri, T rf, int n)
  {
    lower_bound       = ri;
    upper_bound       = rf;
    double ratio      = rf / ri;
    double log_ratio  = std::log(ratio);
    double dlog_ratio = log_ratio / static_cast<double>(n - 1);
    OneOverLogDelta   = 1.0 / dlog_ratio;
    LogDelta          = dlog_ratio;
    r_values.resize(n);
    for (int i = 0; i < n; i++)
      r_values[i] = ri * std::exp(i * dlog_ratio);
  }

  PRAGMA_OFFLOAD("omp declare target")
  static inline T getCL(T r, int& loc, double one_over_log_delta, T lower_bound, double log_delta)
  {
    loc = static_cast<int>(std::log(r / lower_bound) * one_over_log_delta);
    return r - lower_bound * std::exp(loc * log_delta);
  }
  PRAGMA_OFFLOAD("omp end declare target")

  inline int locate(T r) const
  {
    int loc = static_cast<int>(std::log(r / lower_bound) * OneOverLogDelta);
    if (loc >= r_values.size())
      throw std::domain_error("r value " + std::to_string(r) + ">=" + std::to_string(upper_bound) + "\n");
    return loc;
  }

  inline T operator()(int i)
  {
    //return static_cast<T>(lower_bound*std::exp(i*LogDelta));
    return r_values[i];
  }

  //CHECK MIXED PRECISION SENSITIVITY
  inline T getCLForQuintic(T r, int& loc) const
  {
    loc = locate(r);
    return r - r_values[loc];
  }

  inline int size() { return r_values.size(); }
};

/** multivalue implementation for OneDimQuintic
   *  Real valued only
   *  calling any eval method with r >= r_max will throw an exception
   */
template<typename T>
class MultiQuinticSpline1D
{
public:
  using RealType       = T;
  using GridType       = OneDimGridBase<T>;
  using CoeffType      = Matrix<T, OffloadPinnedAllocator<T>>;
  using OffloadArray2D = Array<T, 2, OffloadPinnedAllocator<T>>;
  using OffloadArray3D = Array<T, 3, OffloadPinnedAllocator<T>>;
  using OffloadArray4D = Array<T, 4, OffloadPinnedAllocator<T>>;

private:
  ///number of splines
  size_t num_splines_;
  ///order of spline
  size_t spline_order;
  ///maximum radius
  T r_max;

  ///will be the real grid
  LogGridLight<T> myGrid;

  ///coeffs[6*spline_points][num_splines+padding]
  const std::shared_ptr<CoeffType> coeffs_ptr_;
  const std::shared_ptr<Vector<T, OffloadPinnedAllocator<T>>> first_deriv_ptr_;

  CoeffType& coeffs_;
  Vector<T, OffloadPinnedAllocator<T>>& first_deriv_;

public:
  MultiQuinticSpline1D()
      : coeffs_ptr_(std::make_shared<CoeffType>()),
        first_deriv_ptr_(std::make_shared<Vector<T, OffloadPinnedAllocator<T>>>()),
        coeffs_(*coeffs_ptr_),
        first_deriv_(*first_deriv_ptr_)
  {}

  MultiQuinticSpline1D(const MultiQuinticSpline1D& in) = default;

  inline T rmax() const { return myGrid.upper_bound; }

  inline void evaluate(T r, T* restrict u) const
  {
    if (r < myGrid.lower_bound)
    {
      const T dr          = r - myGrid.lower_bound;
      const T* restrict a = coeffs_[0];
      for (size_t i = 0; i < num_splines_; ++i)
        u[i] = a[i] + first_deriv_[i] * dr;
    }
    else
    {
      int loc;
      const auto cL       = myGrid.getCLForQuintic(r, loc);
      const size_t offset = loc * 6;
      //coeffs is an OhmmsMatrix and [] is a row access operator
      //returning a pointer to 'row' which is normal type pointer []
      const T* restrict a = coeffs_[offset + 0];
      const T* restrict b = coeffs_[offset + 1];
      const T* restrict c = coeffs_[offset + 2];
      const T* restrict d = coeffs_[offset + 3];
      const T* restrict e = coeffs_[offset + 4];
      const T* restrict f = coeffs_[offset + 5];
      for (size_t i = 0; i < num_splines_; ++i)
        u[i] = a[i] + cL * (b[i] + cL * (c[i] + cL * (d[i] + cL * (e[i] + cL * f[i]))));
    }
  }

  /**
   * @brief evaluate MultiQuinticSpline1D for multiple electrons and multiple pbc images
   * 
   * @param [in] r electron distances [Nelec, Npbc]
   * @param [out] u value of all splines at all electron distances [Nelec, Npbc, Nsplines]
   * @param Rmax spline will evaluate to zero for any distance greater than or equal to Rmax
  */
inline void batched_evaluate(const OffloadArray2D& r,
                                   OffloadArray3D& u,
                                   T Rmax) const
{
  // 1) Basic sizes
  const size_t nElec = r.size(0);
  const size_t Nxyz  = r.size(1);
  assert(nElec == u.size(0));
  assert(Nxyz == u.size(1));

  const size_t nRnl  = u.size(2);       // number of splines
  const size_t nR    = nElec * Nxyz;    // total [electron,image] positions
  const size_t nsplines = num_splines_;
  assert(nRnl == nsplines);

  // 2) Obtain device pointers
  const T* r_ptr        = r.device_data();
  T*       u_ptr        = u.device_data();
  const T* coeff_ptr    = coeffs_.device_data();
  const T* fderiv_ptr   = first_deriv_.device_data();

  // 3) Using stack variables to trigger firstprivate transfer  
  const double one_over_log_delta = myGrid.OneOverLogDelta;
  const T      lower_bound        = myGrid.lower_bound;
  const T      log_delta          = myGrid.LogDelta;
  const size_t nCols              = coeffs_.cols();

  // 4) Flatten loops: total threads = (nR * nsplines).
  //    Each thread will handle exactly one (ir, i) pair:
  const size_t total_work = nR * nsplines;

    PRAGMA_OFFLOAD("omp target teams distribute parallel for \
    is_device_ptr(r_ptr, u_ptr, coeff_ptr, fderiv_ptr)")
  for (size_t idx = 0; idx < total_work; idx++)
  {
    // Decode which electron/image index (ir) and which spline (i)
    const size_t ir = idx / nsplines;
    const size_t i  = idx % nsplines;

    // Current radius
    const T rad = r_ptr[ir];

    // Output index in 'u' for this (ir, i) pair
    const size_t outIndex = ir * nRnl + i;

    // Evaluate depending on which region rad is in
    if (rad >= Rmax)
    {
      // Region 1: Above Rmax -> value = 0
      u_ptr[outIndex] = T(0);
    }
    else if (rad < lower_bound)
    {
      // Region 2: Below lower_bound -> linear extrapolation
      const T dr = rad - lower_bound;
      u_ptr[outIndex] = coeff_ptr[i] + fderiv_ptr[i] * dr;
    }
    else
    {
      // Region 3: Main log-spline region
      int loc = 0;
      const T cL = LogGridLight<T>::getCL(rad, loc,
                                          one_over_log_delta,
                                          lower_bound,
                                          log_delta);

      // Each segment has 6 polynomials in the coefficient table
      const size_t offset = loc * 6;
      const T* a = coeff_ptr + nCols * (offset + 0);
      const T* b = coeff_ptr + nCols * (offset + 1);
      const T* c = coeff_ptr + nCols * (offset + 2);
      const T* d = coeff_ptr + nCols * (offset + 3);
      const T* e = coeff_ptr + nCols * (offset + 4);
      const T* f = coeff_ptr + nCols * (offset + 5);

      // Horner-like polynomial evaluation
      u_ptr[outIndex] =
        a[i] + cL * ( b[i] + cL * ( c[i] + cL * ( d[i] +
                   cL * ( e[i] + cL * f[i] ) ) ) );
    }
  }
}

  /**
   * @brief evaluate value, first deriv, second deriv of MultiQuinticSpline1D for multiple electrons and multiple pbc images
   * 
   * r is assumed to be up-to-date on the device when entering this function, and
   * vgl will be up to date on the device when exiting this function
   * 
   * @param [in] r electron distances [Nelec, Npbc]
   * @param [out] vgl val/d1/d2 of all splines at all electron distances [3(val,d1,d2), Nelec, Npbc, Nsplines]
   * @param Rmax spline and derivatives will evaluate to zero for any distance greater than or equal to Rmax
  */
inline void batched_evaluateVGL(const OffloadArray2D& r, OffloadArray4D& vgl, T Rmax) const
{
  const size_t nElec = r.size(0);
  const size_t Nxyz  = r.size(1);
  // vgl shape [3, nElec, Nxyz, nSplines].
  //  "nRnl" is the same as the number of spline
  const size_t nRnl  = vgl.size(3);
  const size_t nR    = nElec * Nxyz; // total positions to evaluate (ir range)

  // Assumption: nRnl == num_splines_, assert(nRnl == num_splines_) :
  const size_t nsplines = num_splines_;
  assert(nRnl == nsplines);

  const T*   r_ptr    = r.device_data();
  T*         u_ptr    = vgl.device_data_at(0, 0, 0, 0);
  T*         du_ptr   = vgl.device_data_at(1, 0, 0, 0);
  T*         d2u_ptr  = vgl.device_data_at(2, 0, 0, 0);
  const T*   coeff_ptr= coeffs_.device_data();
  const T*   fderiv_ptr = first_deriv_.device_data();

  //local copies of frequently used scalars 
  const double one_over_log_delta = myGrid.OneOverLogDelta;
  const T lower_bound             = myGrid.lower_bound;
  const T dlog_ratio              = myGrid.LogDelta;
  const size_t nCols              = coeffs_.cols();

  constexpr T ctwo(2), cthree(3), cfour(4), cfive(5), csix(6), c12(12), c20(20);

  // Flatten loops: total threads = nR * nsplines
  const size_t total_work = nR * nsplines;

    PRAGMA_OFFLOAD("omp target teams distribute parallel for \
    is_device_ptr(r_ptr, u_ptr, du_ptr, d2u_ptr, coeff_ptr, fderiv_ptr)")
  for (size_t idx = 0; idx < total_work; idx++)
  {
    // Decode 2D index from flattened idx
    const size_t ir  = idx / nsplines;  // which electron/image pair
    const size_t i   = idx % nsplines;  // which spline

    // Current radius
    T rad = r_ptr[ir];

    // Prepare local results
    T val, der, d2val;

    if (rad >= Rmax)
    {
      // Past cutoff -> 0
      val   = T(0);
      der   = T(0);
      d2val = T(0);
    }
    else if(rad < lower_bound)
    {

      // Smallest grid point -> linear extrapolation
      const T dr = rad - lower_bound;
      val   = coeff_ptr[i] + fderiv_ptr[i] * dr;
      der   = fderiv_ptr[i];
      d2val = T(0);
    }
    else
    {
      // Log-spline region
      int loc = 0;
      const T cL = LogGridLight<T>::getCL(rad, loc, one_over_log_delta, lower_bound, dlog_ratio);

      // Each spline segment has 6 rows in coeffs_
      // offset picks which polynomial's (a,b,c,d,e,f) row we use
      const size_t offset = loc * 6;

      const T* a = coeff_ptr + nCols * (offset + 0);
      const T* b = coeff_ptr + nCols * (offset + 1);
      const T* c = coeff_ptr + nCols * (offset + 2);
      const T* d = coeff_ptr + nCols * (offset + 3);
      const T* e = coeff_ptr + nCols * (offset + 4);
      const T* f = coeff_ptr + nCols * (offset + 5);

      // Horner-like polynomial expansions
      val = a[i] + cL * (b[i] + cL * (c[i] + cL * (d[i] + cL * (e[i] + cL * f[i]))));

      der = b[i] + cL * (ctwo * c[i] +
                        cL * (cthree * d[i] +
                              cL * (cfour * e[i] + cfive * f[i] * cL)));

      d2val = ctwo * c[i] + cL * (csix * d[i] +
                                 cL * (c12 * e[i] + c20 * f[i] * cL));
    }

    // Store results into the vgl array
    const size_t outIndex = ir * nRnl + i; // same as ir*nsplines+i
    u_ptr[outIndex]   = val;
    du_ptr[outIndex]  = der;
    d2u_ptr[outIndex] = d2val;
  }
}
  inline void evaluate(T r, T* restrict u, T* restrict du, T* restrict d2u) const
  {
    if (r < myGrid.lower_bound)
    {
      const T dr          = r - myGrid.lower_bound;
      const T* restrict a = coeffs_[0];
      for (size_t i = 0; i < num_splines_; ++i)
      {
        u[i]   = a[i] + first_deriv_[i] * dr;
        du[i]  = first_deriv_[i];
        d2u[i] = 0.0;
      }
    }
    else
    {
      int loc;
      const auto cL       = myGrid.getCLForQuintic(r, loc);
      const size_t offset = loc * 6;

      constexpr T ctwo(2);
      constexpr T cthree(3);
      constexpr T cfour(4);
      constexpr T cfive(5);
      constexpr T csix(6);
      constexpr T c12(12);
      constexpr T c20(20);

      const T* restrict a = coeffs_[offset + 0];
      const T* restrict b = coeffs_[offset + 1];
      const T* restrict c = coeffs_[offset + 2];
      const T* restrict d = coeffs_[offset + 3];
      const T* restrict e = coeffs_[offset + 4];
      const T* restrict f = coeffs_[offset + 5];

      for (size_t i = 0; i < num_splines_; ++i)
      {
        u[i]   = a[i] + cL * (b[i] + cL * (c[i] + cL * (d[i] + cL * (e[i] + cL * f[i]))));
        du[i]  = b[i] + cL * (ctwo * c[i] + cL * (cthree * d[i] + cL * (cfour * e[i] + cL * f[i] * cfive)));
        d2u[i] = ctwo * c[i] + cL * (csix * d[i] + cL * (c12 * e[i] + cL * f[i] * c20));
      }
    }
  }

  /** compute upto 3rd derivatives */
  inline void evaluate(T r, T* restrict u, T* restrict du, T* restrict d2u, T* restrict d3u) const
  {
    if (r < myGrid.lower_bound)
    {
      const T dr          = r - myGrid.lower_bound;
      const T* restrict a = coeffs_[0];
      for (size_t i = 0; i < num_splines_; ++i)
      {
        u[i]   = a[i] + first_deriv_[i] * dr;
        du[i]  = first_deriv_[i];
        d2u[i] = 0.0;
        d3u[i] = 0.0;
      }
    }
    else
    {
      int loc;
      const auto cL       = myGrid.getCLForQuintic(r, loc);
      const size_t offset = loc * 6;

      constexpr T ctwo(2);
      constexpr T cthree(3);
      constexpr T cfour(4);
      constexpr T cfive(5);
      constexpr T csix(6);
      constexpr T c12(12);
      constexpr T c20(20);
      constexpr T c24(24);
      constexpr T c60(60);

      const T* restrict a = coeffs_[offset + 0];
      const T* restrict b = coeffs_[offset + 1];
      const T* restrict c = coeffs_[offset + 2];
      const T* restrict d = coeffs_[offset + 3];
      const T* restrict e = coeffs_[offset + 4];
      const T* restrict f = coeffs_[offset + 5];

      for (size_t i = 0; i < num_splines_; ++i)
      {
        u[i]   = a[i] + cL * (b[i] + cL * (c[i] + cL * (d[i] + cL * (e[i] + cL * f[i]))));
        du[i]  = b[i] + cL * (ctwo * c[i] + cL * (cthree * d[i] + cL * (cfour * e[i] + cL * f[i] * cfive)));
        d2u[i] = ctwo * c[i] + cL * (csix * d[i] + cL * (c12 * e[i] + cL * f[i] * c20));
        d3u[i] = csix * d[i] + cL * (c24 * e[i] + cL * f[i] * c60);
      }
    }
  }

  /** initialize grid and container 
   * @param ri minimum  grid point
   * @param rf maximum grid point
   * @param npts number of grid points
   * @param n number of splines
   * @param oreder 5=quintic and 3=cubic
   */
  template<typename GT>
  void initialize(GT& agrid, int norbs, int order = 5)
  {
    myGrid.set(agrid.rmin(), agrid.rmax(), agrid.size());
    r_max = myGrid.upper_bound;
    if (coeffs_.size())
      throw std::runtime_error("MultiQuinticSpline1D::initialize cannot reinitialize coeffs.");

    spline_order = order;
    num_splines_ = norbs;
    coeffs_.resize((order + 1) * agrid.size(), getAlignedSize<T>(norbs));
    first_deriv_.resize(num_splines_);
  }

  template<typename T1>
  void add_spline(int ispline, OneDimQuinticSpline<T1>& in)
  {
    first_deriv_[ispline] = in.first_deriv;
    //if(spline_order==QUINTIC)
    {
      const T1* restrict A    = in.m_Y.data();
      const T1* restrict B    = in.B.data();
      const T1* restrict C    = in.m_Y2.data();
      const T1* restrict D    = in.D.data();
      const T1* restrict E    = in.E.data();
      const T1* restrict F    = in.F.data();
      T* restrict out         = coeffs_.data();
      const size_t ncols      = coeffs_.cols();
      const size_t num_points = in.size();
      for (size_t i = 0; i < num_points; ++i)
      {
        out[(i * 6 + 0) * ncols + ispline] = static_cast<T>(A[i]);
        out[(i * 6 + 1) * ncols + ispline] = static_cast<T>(B[i]);
        out[(i * 6 + 2) * ncols + ispline] = static_cast<T>(C[i]);
        out[(i * 6 + 3) * ncols + ispline] = static_cast<T>(D[i]);
        out[(i * 6 + 4) * ncols + ispline] = static_cast<T>(E[i]);
        out[(i * 6 + 5) * ncols + ispline] = static_cast<T>(F[i]);
      }
    }
  }

  /// finalize the construction
  void finalize()
  {
    first_deriv_.updateTo();
    coeffs_.updateTo();
  }

  int getNumSplines() const { return num_splines_; }
  void setNumSplines(int num_splines) { num_splines_ = num_splines; }
};

extern template class MultiQuinticSpline1D<float>;
extern template class MultiQuinticSpline1D<double>;
} // namespace qmcplusplus
#endif
