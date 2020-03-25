////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineData.hpp
 *
 * header file for shared service
 */
#ifndef SPLINE2_MULTIEINSPLINE_DATA_HPP
#define SPLINE2_MULTIEINSPLINE_DATA_HPP

namespace spline2
{
/** class for cublic spline parameters
 *
 * static constant data and static functions are defined here
 */
template<typename T>
struct MultiBsplineData
{
  // clang-format off
  /// spline interpolation parameters, statically initialized
  static constexpr T A00 = -1.0/6.0, A01 =  3.0/6.0, A02 = -3.0/6.0, A03 = 1.0/6.0;
  static constexpr T A10 =  3.0/6.0, A11 = -6.0/6.0, A12 =  0.0/6.0, A13 = 4.0/6.0;
  static constexpr T A20 = -3.0/6.0, A21 =  3.0/6.0, A22 =  3.0/6.0, A23 = 1.0/6.0;
  static constexpr T A30 =  1.0/6.0, A31 =  0.0/6.0, A32 =  0.0/6.0, A33 = 0.0/6.0;
  static constexpr T dA01 = -0.5, dA02 =  1.0, dA03 = -0.5;
  static constexpr T dA11 =  1.5, dA12 = -2.0, dA13 =  0.0;
  static constexpr T dA21 = -1.5, dA22 =  1.0, dA23 =  0.5;
  static constexpr T dA31 =  0.5, dA32 =  0.0, dA33 =  0.0;
  static constexpr T d2A02 = -1.0, d2A03 =  1.0;
  static constexpr T d2A12 =  3.0, d2A13 = -2.0;
  static constexpr T d2A22 = -3.0, d2A23 =  1.0;
  static constexpr T d2A32 =  1.0, d2A33 =  0.0;
  static constexpr T d3A03 = -1.0;
  static constexpr T d3A13 =  3.0;
  static constexpr T d3A23 = -3.0;
  static constexpr T d3A33 =  1.0;
  // clang-format on

  /** compute interpolation prefactors
   * @param tx fractional coordinate with respect to the grid length in range [0,1)
   * @param a[4] prefactor for the four consecutive grid points
   */
  inline static void compute_prefactors(T a[4], T tx)
  {
    a[0] = ((A00 * tx + A01) * tx + A02) * tx + A03;
    a[1] = ((A10 * tx + A11) * tx + A12) * tx + A13;
    a[2] = ((A20 * tx + A21) * tx + A22) * tx + A23;
    a[3] = ((A30 * tx + A31) * tx + A32) * tx + A33;
  }

  /** compute interpolation prefactors up to the second order
   * @param tx fractional coordinate with respect to the grid length in range [0,1)
   * @param a[4] prefactor for the four consecutive grid points
   * @param da[4] first order prefactor for the four consecutive grid points
   * @param d2a[4] second order prefactor for the four consecutive grid points
   */
  inline static void compute_prefactors(T a[4], T da[4], T d2a[4], T tx)
  {
    a[0]   = ((A00 * tx + A01) * tx + A02) * tx + A03;
    a[1]   = ((A10 * tx + A11) * tx + A12) * tx + A13;
    a[2]   = ((A20 * tx + A21) * tx + A22) * tx + A23;
    a[3]   = ((A30 * tx + A31) * tx + A32) * tx + A33;
    da[0]  = (dA01 * tx + dA02) * tx + dA03;
    da[1]  = (dA11 * tx + dA12) * tx + dA13;
    da[2]  = (dA21 * tx + dA22) * tx + dA23;
    da[3]  = (dA31 * tx + dA32) * tx + dA33;
    d2a[0] = d2A02 * tx + d2A03;
    d2a[1] = d2A12 * tx + d2A13;
    d2a[2] = d2A22 * tx + d2A23;
    d2a[3] = d2A32 * tx + d2A33;
  }

  /** compute interpolation prefactors up to the third order
   * @param tx fractional coordinate with respect to the grid length in range [0,1)
   * @param a[4] prefactor for the four consecutive grid points
   * @param da[4] first order prefactor for the four consecutive grid points
   * @param d2a[4] second order prefactor for the four consecutive grid points
   * @param d3a[4] third order prefactor for the four consecutive grid points
   */
  inline static void compute_prefactors(T a[4], T da[4], T d2a[4], T d3a[4], T tx)
  {
    a[0]   = ((A00 * tx + A01) * tx + A02) * tx + A03;
    a[1]   = ((A10 * tx + A11) * tx + A12) * tx + A13;
    a[2]   = ((A20 * tx + A21) * tx + A22) * tx + A23;
    a[3]   = ((A30 * tx + A31) * tx + A32) * tx + A33;
    da[0]  = (dA01 * tx + dA02) * tx + dA03;
    da[1]  = (dA11 * tx + dA12) * tx + dA13;
    da[2]  = (dA21 * tx + dA22) * tx + dA23;
    da[3]  = (dA31 * tx + dA32) * tx + dA33;
    d2a[0] = d2A02 * tx + d2A03;
    d2a[1] = d2A12 * tx + d2A13;
    d2a[2] = d2A22 * tx + d2A23;
    d2a[3] = d2A32 * tx + d2A33;
    d3a[0] = d3A03;
    d3a[1] = d3A13;
    d3a[2] = d3A23;
    d3a[3] = d3A33;
  }
};

} // namespace spline2

#endif
