//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_BSPLINE_TRAITS_H
#define QMCPLUSPLUS_BSPLINE_TRAITS_H

#include <complex>
#include "einspline/multi_bspline.h"
namespace qmcplusplus
{
/** determine if EngT (e.g., einspline engine) handles real data or complex data
   *
   * Default is true and complex cases are specialized
   */
template<typename EngT>
struct is_real_bspline
{
  static const bool value = true;
};

///specialization for multi_UBspline_3d_z
template<>
struct is_real_bspline<multi_UBspline_3d_z>
{
  static const bool value = false;
};

///specialization for multi_UBspline_3d_z
template<>
struct is_real_bspline<multi_UBspline_3d_c>
{
  static const bool value = false;
};
//
/** struct to check if two types are the same
   */
template<typename T1, typename T2>
struct type_check
{
  static const bool value = false;
};

template<typename T1>
struct type_check<T1, T1>
{
  static const bool value = true;
};

template<typename T1, typename T2>
struct is_same_precision
{
  static const bool value = false;
};
template<>
struct is_same_precision<double, double>
{
  static const bool value = true;
};
template<>
struct is_same_precision<float, float>
{
  static const bool value = true;
};

/** dummy traits class for bspline engine
   *
   * Should fail to instantiate invalid engines if the trait class is not implemented
   * The specialization provides
   * - DIM, physical dimension
   * - real_type real data type
   * - value_type value data type
   * - spline_type einspline engine type
   */
template<typename EngT>
struct bspline_engine_traits
{};

#if OHMMS_DIM == 3
/** specialization with multi_UBspline_3d_d
   */
template<>
struct bspline_engine_traits<multi_UBspline_3d_d>
{
  enum
  {
    DIM = 3
  };
  using SplineType       = multi_UBspline_3d_d;
  using SingleSplineType = UBspline_3d_d;
  using BCType           = BCtype_d;
  using real_type        = double;
  using value_type       = double;
};

///specialization with multi_UBspline_3d_z
template<>
struct bspline_engine_traits<multi_UBspline_3d_z>
{
  enum
  {
    DIM = 3
  };
  using SplineType         = multi_UBspline_3d_z;
  using SingleSplineType   = UBspline_3d_z;
  using BCType             = BCtype_z;
  using value_type         = std::complex<double>;
  using single_spline_type = UBspline_3d_z;
};

/** specialization with multi_UBspline_3d_s
   */
template<>
struct bspline_engine_traits<multi_UBspline_3d_s>
{
  enum
  {
    DIM = 3
  };
  using SplineType       = multi_UBspline_3d_s;
  using SingleSplineType = UBspline_3d_s;
  using BCType           = BCtype_s;
  using real_type        = float;
  using value_type       = float;
};

///specialization with multi_UBspline_3d_c
template<>
struct bspline_engine_traits<multi_UBspline_3d_c>
{
  enum
  {
    DIM = 3
  };
  using SplineType       = multi_UBspline_3d_c;
  using SingleSplineType = UBspline_3d_c;
  using BCType           = BCtype_c;
  using real_type        = float;
  using value_type       = std::complex<float>;
};

#else
/** specialization with multi_UBspline_3d_d
   */
template<>
struct bspline_engine_traits<multi_UBspline_2d_d>
{
  enum
  {
    DIM = 2
  };
  using real_type  = double;
  using value_type = double;
};

///specialization with multi_UBspline_3d_z
template<>
struct bspline_engine_traits<multi_UBspline_2d_z>
{
  enum
  {
    DIM = 2
  };
  using real_type  = double;
  using value_type = std::complex<double>;
};

/** specialization with multi_UBspline_3d_d
   */
template<>
struct bspline_engine_traits<multi_UBspline_2d_s>
{
  enum
  {
    DIM = 2
  };
  using real_type  = float;
  using value_type = float;
};

///specialization with multi_UBspline_3d_z
template<>
struct bspline_engine_traits<multi_UBspline_2d_c>
{
  enum
  {
    DIM = 2
  };
  using real_type  = float;
  using value_type = std::complex<float>;
};
#endif

/** specialization with multi_UBspline_3d_d
   */
template<>
struct bspline_engine_traits<multi_UBspline_1d_d>
{
  enum
  {
    DIM = 1
  };
  using SplineType       = multi_UBspline_1d_d;
  using SingleSplineType = UBspline_1d_d;
  using BCType           = BCtype_d;
  using real_type        = double;
  using value_type       = double;
};

/** specialization with multi_UBspline_3d_s
   */
template<>
struct bspline_engine_traits<multi_UBspline_1d_s>
{
  enum
  {
    DIM = 1
  };
  using SplineType       = multi_UBspline_1d_s;
  using SingleSplineType = UBspline_1d_s;
  using BCType           = BCtype_s;
  using real_type        = float;
  using value_type       = float;
};

} // namespace qmcplusplus
#endif
