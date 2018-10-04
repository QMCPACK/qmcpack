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

#include <cmath>
#include <algorithm>

namespace spline2
{

template<typename T, typename TRESIDUAL>
inline void getSplineBound(T x, TRESIDUAL& dx, int& ind, int ng)
{
  T ipart;
  dx=std::modf(x,&ipart);
  ind = std::min(std::max(int(0),static_cast<int>(ipart)),ng);
}

template<typename T>
struct MultiBsplineData
{
  static const T A44[16];
  static const T dA44[16];
  static const T d2A44[16];
  static const T d3A44[16];

  inline static void compute_prefactors(T a[4], T tx)
  {
    a[0] = ((A44[0] * tx + A44[1]) * tx + A44[2]) * tx + A44[3];
    a[1] = ((A44[4] * tx + A44[5]) * tx + A44[6]) * tx + A44[7];
    a[2] = ((A44[8] * tx + A44[9]) * tx + A44[10]) * tx + A44[11];
    a[3] = ((A44[12] * tx + A44[13]) * tx + A44[14]) * tx + A44[15];
  }

  inline static void compute_prefactors(T a[4], T da[4], T d2a[4], T tx)
  {
    a[0]   = ((A44[0] * tx + A44[1]) * tx + A44[2]) * tx + A44[3];
    a[1]   = ((A44[4] * tx + A44[5]) * tx + A44[6]) * tx + A44[7];
    a[2]   = ((A44[8] * tx + A44[9]) * tx + A44[10]) * tx + A44[11];
    a[3]   = ((A44[12] * tx + A44[13]) * tx + A44[14]) * tx + A44[15];
    da[0]  = ((dA44[0] * tx + dA44[1]) * tx + dA44[2]) * tx + dA44[3];
    da[1]  = ((dA44[4] * tx + dA44[5]) * tx + dA44[6]) * tx + dA44[7];
    da[2]  = ((dA44[8] * tx + dA44[9]) * tx + dA44[10]) * tx + dA44[11];
    da[3]  = ((dA44[12] * tx + dA44[13]) * tx + dA44[14]) * tx + dA44[15];
    d2a[0] = ((d2A44[0] * tx + d2A44[1]) * tx + d2A44[2]) * tx + d2A44[3];
    d2a[1] = ((d2A44[4] * tx + d2A44[5]) * tx + d2A44[6]) * tx + d2A44[7];
    d2a[2] = ((d2A44[8] * tx + d2A44[9]) * tx + d2A44[10]) * tx + d2A44[11];
    d2a[3] = ((d2A44[12] * tx + d2A44[13]) * tx + d2A44[14]) * tx + d2A44[15];
  }
};

} // namespace spline2

#endif
