////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineEvalHelper.hpp
 *
 * header file for helper functions
 */
#ifndef SPLINE2_MULTIEINSPLINE_EVAL_HELPER_HPP
#define SPLINE2_MULTIEINSPLINE_EVAL_HELPER_HPP

#include <cmath>
#include <algorithm>

namespace spline2
{

/** break x into the integer part and residual part and apply bounds
 * @param x input coordinate
 * @param dx fractional part
 * @param ind integer part
 * @param nmax upper bound of the integer part
 *
 * x in the range of [0, nmax+1) will be split correctly.
 * x < 0, ind = 0, dx = 0
 * x >= nmax+1, ind = nmax, dx = 1 - epsilon
 *
 * Attention: nmax is not the number grid points but the maximum allowed grid index
 * For example, ng is the number of grid point.
 * the actual grid points indices are 0, 1, ..., ng - 1.
 * In a periodic/anti periodic spline, set nmax = ng - 1
 * In a natural boundary spline, set nmax = ng - 2
 * because the end point should be excluded and the last grid point has an index ng - 2.
 */
template<typename T, typename TRESIDUAL>
inline void getSplineBound(T x, TRESIDUAL& dx, int& ind, int nmax)
{
  // lower bound
  if (x < 0)
  {
    ind = 0;
    dx  = T(0);
  }
  else
  {
    T ipart;
    dx  = std::modf(x, &ipart);
    ind = static_cast<int>(ipart);
    // upper bound
    if (ind > nmax)
    {
      ind = nmax;
      dx  = T(1) - std::numeric_limits<T>::epsilon();
    }
  }
}

} // namespace spline2

#endif
