//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LINEAR_GRID_CUBIC_SPLINE_H
#define QMCPLUSPLUS_LINEAR_GRID_CUBIC_SPLINE_H

#include "OneDimCubicSpline.h"
#include "OneDimCubicSpline.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

namespace qmcplusplus
{
/** combined OneDimCubicSpline and LinearGrid
 * OneDimCubicSpline contains OneDimGridBase pointer and calls its virtual function members. This doesn't work well on a GPU.
 * Since the use case is OneDimCubicSpline with LinearGrid. We fuse both classes and avoid any virtual functions.
 * There are two splint functions. The one with one paramemter r is intended for testing or being called on the CPU.
 * The static one with many parameters is intended to be used(inlined) inside a GPU kernel.
 */
template<class T>
class OneDimCubicSplineLinearGrid
{
public:
  OneDimCubicSplineLinearGrid(const OneDimCubicSpline<T>& cublis_spliner)
  {
    auto& grid = dynamic_cast<const LinearGrid<T>&>(cublis_spliner.grid());
    r_min_     = grid.rmin();
    r_max_     = grid.rmax();
    delta_inv_ = grid.DeltaInv;

    const size_t num_points = grid.size();
    X_.resize(num_points);
    for (size_t i = 0; i < num_points; i++)
      X_[i] = *(grid.data() + i);

    first_deriv_ = cublis_spliner.first_deriv;
    const_value_ = cublis_spliner.ConstValue;
    m_Y_.resize(num_points);
    m_Y2_.resize(num_points);
    for (size_t i = 0; i < num_points; i++)
    {
      m_Y_[i]  = cublis_spliner.m_Y[i];
      m_Y2_[i] = cublis_spliner.m_Y2[i];
    }
    X_.updateTo();
    m_Y_.updateTo();
    m_Y2_.updateTo();
  }

  /** compute the function value at r
   */
  T splint(T r) const
  {
    return splint(r_min_, r_max_, X_.data(), delta_inv_, m_Y_.data(), m_Y2_.data(), first_deriv_, const_value_, r);
  }

  /** compute the function value at r.
   * Need to pass in all the parameters.
   */
  static T splint(T r_min,
                  T r_max,
                  const T* X,
                  T delta_inv,
                  const T* m_Y,
                  const T* m_Y2,
                  T first_deriv,
                  T const_value,
                  T r)
  {
    if (r < r_min)
    {
      return m_Y[0] + first_deriv * (r - r_min);
    }
    else if (r >= r_max)
    {
      return const_value;
    }

    const size_t loc = std::floor((r - r_min) * delta_inv);
    const T dist     = r - X[loc];
    const T delta    = X[loc + 1] - X[loc];
    CubicSplineEvaluator<T> eval(dist, delta);
    return eval.cubicInterpolate(m_Y[loc], m_Y[loc + 1], m_Y2[loc], m_Y2[loc + 1]);
  }

  const auto& get_m_Y() const { return m_Y_; }
  const auto& get_m_Y2() const { return m_Y2_; }
  T get_first_deriv() const { return first_deriv_; }
  T get_const_value() const { return const_value_; }

  T get_r_min() const { return r_min_; }
  T get_r_max() const { return r_max_; }
  const auto& get_X() const { return X_; }
  double get_delta_inv() const { return delta_inv_; }

private:
  // spline related
  /// data for the function on the grid
  Vector<T, OffloadAllocator<T>> m_Y_;
  /// data for the function on the grid
  Vector<T, OffloadAllocator<T>> m_Y2_;
  /// first derivative for handling r < r_min_
  T first_deriv_;
  /// const value for handling r > r_max_
  T const_value_;

  // grid related info
  /// use spline above r_min_. If below, use first deriv extrapolation
  T r_min_;
  /// use spline below r_min_. If above, use const value
  T r_max_;
  /// the location of grid points
  Vector<T, OffloadAllocator<T>> X_;
  /// 1/grid space
  double delta_inv_;
};

} // namespace qmcplusplus
#endif
