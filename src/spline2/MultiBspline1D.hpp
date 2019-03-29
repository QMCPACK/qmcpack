//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBspline.hpp
 *
 * define classes MultiBspline1D
 */

#ifndef QMCPLUSPLUS_MULTIEINSPLINE_1D_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_1D_HPP

namespace qmcplusplus
{
/** container class to hold a 1D multi spline structure
   * @tparam T the precision of splines
   */
template<typename T>
struct MultiBspline1D
{
  ///define the einsplie object type
  using SplineType = typename bspline_traits<T, 1>::SplineType;
  ///define the real type
  using real_type = typename bspline_traits<T, 1>::real_type;
  ///actual einspline multi-bspline object
  SplineType spline_m;

  MultiBspline1D()
  {
    spline_m.coefs       = nullptr;
    spline_m.num_splines = 0;
    spline_m.coefs_size  = 0;
  }

  /** create the einspline as used in the builder
       */
  template<typename GT, typename BCT>
  void create(GT& grid, BCT& bc, int num_splines)
  {
    if (getAlignedSize<T>(num_splines) != num_splines)
      throw std::runtime_error("When creating the data space of MultiBspline1D, num_splines must be padded!\n");
    SplineType* temp_spline;
    temp_spline = einspline::create(temp_spline, grid, bc, num_splines);
    spline_m    = *temp_spline;
    free(temp_spline);
  }

  void flush_zero() const
  {
    if (spline_m.coefs != nullptr)
      std::fill(spline_m.coefs, spline_m.coefs + spline_m.coefs_size, T(0));
  }

  int num_splines() const { return spline_m.num_splines; }

  size_t sizeInByte() const { return (spline_m.coefs == nullptr) ? 0 : spline_m.coefs_size * sizeof(T); }

  /** copy a single spline to the big table
       * @param aSpline UBspline_3d_(d,s)
       * @param int index of aSpline
       * @param offset_ starting index for the case of multiple domains
       * @param base_ number of bases
       */
  template<typename SingleSpline>
  void copy_spline(SingleSpline* aSpline, int i, const int offset_, const int base_)
  {
    einspline::set(&spline_m, i, aSpline, offset_, base_);
  }

  template<typename PT, typename VT>
  inline void evaluate(const PT& r, VT& psi) const
  {
    evaluate_v_impl(r, psi.data());
  }

  template<typename PT, typename VT, typename GT, typename LT>
  inline void evaluate_vgl(const PT& r, VT& psi, GT& grad, LT& lap) const
  {
    evaluate_vgl_impl(r, psi.data(), grad.data(), lap.data());
  }

  //template<typename PT, typename VT, typename GT, typename HT>
  //  inline void evaluate_vgh(const PT& r, VT& psi, GT& grad, HT& hess)
  //  {
  //    evaluate_vgh_impl(r,psi.data(),grad.data(),hess.data());
  //  }

  /// compute values only.
  void evaluate_v_impl(T r, T* restrict vals) const;
  /// compute VGL.
  void evaluate_vgl_impl(T r, T* restrict vals, T* restrict grads, T* restrict lapl) const;
};

template<typename T>
inline void MultiBspline1D<T>::evaluate_v_impl(T x, T* restrict vals) const
{
  x -= spline_m.x_grid.start;
  T tx;
  int ix;
  spline2::getSplineBound(x * spline_m.x_grid.delta_inv, tx, ix, spline_m.x_grid.num - 2);

  T a[4];
  spline2::MultiBsplineData<T>::compute_prefactors(a, tx);

  const intptr_t xs = spline_m.x_stride;
  const T dxInv     = spline_m.x_grid.delta_inv;

  const T* restrict coefs0 = spline_m.coefs + ((ix)*xs);
  const T* restrict coefs1 = spline_m.coefs + ((ix + 1) * xs);
  const T* restrict coefs2 = spline_m.coefs + ((ix + 2) * xs);
  const T* restrict coefs3 = spline_m.coefs + ((ix + 3) * xs);

#pragma omp simd aligned(vals, coefs0, coefs1, coefs2, coefs3)
  for (int n = 0; n < spline_m.num_splines; n++)
    vals[n] = a[0] * coefs0[n] + a[1] * coefs1[n] + a[2] * coefs2[n] + a[3] * coefs3[n];
}

template<typename T>
inline void MultiBspline1D<T>::evaluate_vgl_impl(T x, T* restrict vals, T* restrict grads, T* restrict lapl) const
{
  x -= spline_m.x_grid.start;
  T tx;
  int ix;
  spline2::getSplineBound(x * spline_m.x_grid.delta_inv, tx, ix, spline_m.x_grid.num - 2);

  T a[4], da[4], d2a[4];
  spline2::MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);

  const intptr_t xs = spline_m.x_stride;
  const T dxInv     = spline_m.x_grid.delta_inv;

  const T* restrict coefs0 = spline_m.coefs + ((ix)*xs);
  const T* restrict coefs1 = spline_m.coefs + ((ix + 1) * xs);
  const T* restrict coefs2 = spline_m.coefs + ((ix + 2) * xs);
  const T* restrict coefs3 = spline_m.coefs + ((ix + 3) * xs);

#pragma omp simd aligned(vals, grads, lapl, coefs0, coefs1, coefs2, coefs3)
  for (int n = 0; n < spline_m.num_splines; n++)
  {
    const T coef_0 = coefs0[n];
    const T coef_1 = coefs1[n];
    const T coef_2 = coefs2[n];
    const T coef_3 = coefs3[n];
    vals[n]        = a[0] * coef_0 + a[1] * coef_1 + a[2] * coef_2 + a[3] * coef_3;
    grads[n]       = (da[0] * coef_0 + da[1] * coef_1 + da[2] * coef_2 + da[3] * coef_3) * dxInv;
    lapl[n]        = (d2a[0] * coef_0 + d2a[1] * coef_1 + d2a[2] * coef_2 + d2a[3] * coef_3) * dxInv * dxInv;
  }
}

} // namespace qmcplusplus
#endif
