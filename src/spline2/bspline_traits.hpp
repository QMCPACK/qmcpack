//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file bspline_traits.h
 *
 * extend spline/bspline_traints by introducing
 * bspline_traits<T,D> with only specialization with 3D
 */
#ifndef QMCPLUSPLUS_BSPLINE_SPLINE2_TRAITS_H
#define QMCPLUSPLUS_BSPLINE_SPLINE2_TRAITS_H


#include "einspline/bspline_base.h"
#include "einspline/bspline_structs.h"
#include "einspline/multi_bspline_structs.h"

namespace qmcplusplus
{
/** trait class to map (datatype,D) to Einspline engine type */
template<typename T, unsigned D>
struct bspline_traits
{};

template<>
struct bspline_traits<float, 3>
{
  using SplineType                       = multi_UBspline_3d_s;
  using SingleSplineType                 = UBspline_3d_s;
  using BCType                           = BCtype_s;
  using real_type                        = float;
  using value_type                       = float;
  static const spline_code spcode        = MULTI_U3D;
  static const spline_code single_spcode = U3D;
  static const type_code tcode           = SINGLE_REAL;
};

template<>
struct bspline_traits<double, 3>
{
  using SplineType                       = multi_UBspline_3d_d;
  using SingleSplineType                 = UBspline_3d_d;
  using BCType                           = BCtype_d;
  using real_type                        = double;
  using value_type                       = double;
  static const spline_code spcode        = MULTI_U3D;
  static const spline_code single_spcode = U3D;
  static const type_code tcode           = DOUBLE_REAL;
};

/** specialization for 1D float */
template<>
struct bspline_traits<float, 1>
{
  using SplineType       = multi_UBspline_1d_s;
  using SingleSplineType = UBspline_1d_s;
  using BCType           = BCtype_s;
  using DataType         = float;
  using real_type        = float;
  using value_type       = float;
};

/** specialization for 1D double */
template<>
struct bspline_traits<double, 1>
{
  using SplineType       = multi_UBspline_1d_d;
  using SingleSplineType = UBspline_1d_d;
  using BCType           = BCtype_d;
  using DataType         = double;
  using real_type        = double;
  using value_type       = double;
};


/** helper class to determine the value_type of einspline objects
   */
template<typename ST>
struct bspline_type
{};

template<>
struct bspline_type<multi_UBspline_3d_s>
{
  using value_type = float;
};

template<>
struct bspline_type<multi_UBspline_3d_d>
{
  using value_type = double;
};

template<>
struct bspline_type<UBspline_3d_s>
{
  using value_type = float;
};

template<>
struct bspline_type<UBspline_3d_d>
{
  using value_type = double;
};
} // namespace qmcplusplus
#endif
