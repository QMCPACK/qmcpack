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
  inline void batched_evaluate(const OffloadArray2D& r, OffloadArray3D& u, T Rmax) const
  {
    const size_t nElec = r.size(0);
    const size_t Nxyz  = r.size(1); // number of PBC images
    assert(nElec == u.size(0));
    assert(Nxyz == u.size(1));
    const size_t nRnl = u.size(2);    // number of splines
    const size_t nR   = nElec * Nxyz; // total number of positions to evaluate

    double one_over_log_delta = myGrid.OneOverLogDelta;
    T lower_bound             = myGrid.lower_bound;
    T log_delta               = myGrid.LogDelta;

    auto* r_ptr = r.data();
    auto* u_ptr = u.data();

    auto* coeff_ptr       = coeffs_.data();
    auto* first_deriv_ptr = first_deriv_.data();
    const size_t nCols    = coeffs_.cols();
    const size_t coefsize = coeffs_.size();
    const int nsplines   = num_splines_;

    PRAGMA_OFFLOAD("omp target teams distribute parallel for \
                    map(to:coeff_ptr[:coefsize], first_deriv_ptr[:nsplines]) \
                    map(to:r_ptr[:nR], u_ptr[:nRnl*nR])")
    for (int ir = 0; ir < nR; ir++)
    {
      if (r_ptr[ir] >= Rmax)
      {
        for (int i = 0; i < nsplines; ++i)
          u_ptr[ir * nRnl + i] = 0.0;
      }
      else if (r_ptr[ir] < lower_bound)
      {
        const T dr = r_ptr[ir] - lower_bound;
        for (int i = 0; i < nsplines; ++i)
          u_ptr[ir * nRnl + i] = coeff_ptr[i] + first_deriv_ptr[i] * dr;
      }
      else
      {
        int loc;
        const auto cL       = LogGridLight<T>::getCL(r_ptr[ir], loc, one_over_log_delta, lower_bound, log_delta);
        const size_t offset = loc * 6;
        const T* restrict a = coeff_ptr + nCols * (offset + 0);
        const T* restrict b = coeff_ptr + nCols * (offset + 1);
        const T* restrict c = coeff_ptr + nCols * (offset + 2);
        const T* restrict d = coeff_ptr + nCols * (offset + 3);
        const T* restrict e = coeff_ptr + nCols * (offset + 4);
        const T* restrict f = coeff_ptr + nCols * (offset + 5);
        for (int i = 0; i < nsplines; ++i)
          u_ptr[ir * nRnl + i] = a[i] + cL * (b[i] + cL * (c[i] + cL * (d[i] + cL * (e[i] + cL * f[i]))));
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
    const size_t Nxyz  = r.size(1); // number of PBC images
    assert(3 == vgl.size(0));
    assert(nElec == vgl.size(1));
    assert(Nxyz == vgl.size(2));
    const size_t nRnl = vgl.size(3);  // number of splines
    const size_t nR   = nElec * Nxyz; // total number of positions to evaluate

    double one_over_log_delta = myGrid.OneOverLogDelta;
    T lower_bound             = myGrid.lower_bound;
    T dlog_ratio              = myGrid.LogDelta;

    auto* r_ptr   = r.data();
    auto* u_ptr   = vgl.data_at(0, 0, 0, 0);
    auto* du_ptr  = vgl.data_at(1, 0, 0, 0);
    auto* d2u_ptr = vgl.data_at(2, 0, 0, 0);

    auto* coeff_ptr       = coeffs_.data();
    auto* first_deriv_ptr = first_deriv_.data();
    const size_t nCols    = coeffs_.cols();
    const size_t coefsize = coeffs_.size();
    const auto nsplines   = num_splines_;

    constexpr T ctwo(2);
    constexpr T cthree(3);
    constexpr T cfour(4);
    constexpr T cfive(5);
    constexpr T csix(6);
    constexpr T c12(12);
    constexpr T c20(20);

    PRAGMA_OFFLOAD("omp target teams distribute parallel for \
                    map(to: first_deriv_ptr[:nsplines], coeff_ptr[:coefsize]) \
                    map(to: r_ptr[:nR], u_ptr[:nRnl*nR], du_ptr[:nRnl*nR], d2u_ptr[:nRnl*nR])")
    for (size_t ir = 0; ir < nR; ir++)
    {
      if (r_ptr[ir] >= Rmax)
      {
        for (size_t i = 0; i < nsplines; ++i)
        {
          u_ptr[ir * nRnl + i]   = 0.0;
          du_ptr[ir * nRnl + i]  = 0.0;
          d2u_ptr[ir * nRnl + i] = 0.0;
        }
      }
      else if (r_ptr[ir] < lower_bound)
      {
        const T dr          = r_ptr[ir] - lower_bound;
        const T* restrict a = coeff_ptr;
        // const T* restrict a = coeffs_[0];
        for (size_t i = 0; i < nsplines; ++i)
        {
          u_ptr[ir * nRnl + i]   = a[i] + first_deriv_ptr[i] * dr;
          du_ptr[ir * nRnl + i]  = first_deriv_ptr[i];
          d2u_ptr[ir * nRnl + i] = 0.0;
        }
      }
      else
      {
        int loc;
        const auto cL = LogGridLight<T>::getCL(r_ptr[ir], loc, one_over_log_delta, lower_bound, dlog_ratio);
        // const auto cL       = myGrid.getCLForQuintic(r_list[ir], loc);
        const size_t offset = loc * 6;
        const T* restrict a = coeff_ptr + nCols * (offset + 0);
        const T* restrict b = coeff_ptr + nCols * (offset + 1);
        const T* restrict c = coeff_ptr + nCols * (offset + 2);
        const T* restrict d = coeff_ptr + nCols * (offset + 3);
        const T* restrict e = coeff_ptr + nCols * (offset + 4);
        const T* restrict f = coeff_ptr + nCols * (offset + 5);
        for (size_t i = 0; i < nsplines; ++i)
        {
          u_ptr[ir * nRnl + i] = a[i] + cL * (b[i] + cL * (c[i] + cL * (d[i] + cL * (e[i] + cL * f[i]))));
          du_ptr[ir * nRnl + i] =
              b[i] + cL * (ctwo * c[i] + cL * (cthree * d[i] + cL * (cfour * e[i] + cL * f[i] * cfive)));
          d2u_ptr[ir * nRnl + i] = ctwo * c[i] + cL * (csix * d[i] + cL * (c12 * e[i] + cL * f[i] * c20));
        }
      }
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
