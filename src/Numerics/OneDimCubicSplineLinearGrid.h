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
 */
template<class T>
struct OneDimCubicSplineLinearGrid
{
  // spline related
  Vector<T, OffloadAllocator<T>> m_Y_;
  Vector<T, OffloadAllocator<T>> m_Y2_;
  T first_deriv_;
  T const_value_;

  // grid related info
  T r_min_;
  T r_max_;
  Vector<T, OffloadAllocator<T>> X_;
  double delta_inv_;

  OneDimCubicSplineLinearGrid(const OneDimCubicSpline<T>& cublis_spliner)
  {
    auto& grid = dynamic_cast<const LinearGrid<T>&>(cublis_spliner.grid());
    r_min_ = grid.rmin();
    r_max_ = grid.rmax();
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
      m_Y_[i] = cublis_spliner.m_Y[i];
      m_Y2_[i] = cublis_spliner.m_Y2[i];
    }
    X_.updateTo();
    m_Y_.updateTo();
    m_Y2_.updateTo();
  }

  T splint(T r) const
  {
    return splint(r_min_, r_max_, X_.data(), delta_inv_, m_Y_.data(), m_Y2_.data(), first_deriv_, const_value_, r);
  }

  static T splint(T r_min, T r_max, const T* X, T delta_inv, const T* m_Y, const T* m_Y2, T first_deriv, T const_value, T r)
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
    const T dist = r - X[loc];
    const T delta = X[loc + 1] - X[loc];
    CubicSplineEvaluator<value_type> eval(dist, delta);
    return eval.cubicInterpolate(m_Y[loc], m_Y[loc + 1], m_Y2[loc], m_Y2[loc + 1]);
  }
};

}
#endif
