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
 * bspline_traits<T,D> with only speciliazation with 3D
 */
#ifndef QMCPLUSPLUS_BSPLINE_SPLINE2_TRAITS_H
#define QMCPLUSPLUS_BSPLINE_SPLINE2_TRAITS_H
#include <einspline/multi_bspline.h>
namespace qmcplusplus
{
  /** trait class to map (datatype,D) to Einspline engine type */
  template<typename T, unsigned D> struct bspline_traits{};

  template<>
    struct bspline_traits<float,3>
    {
      typedef multi_UBspline_3d_s SplineType;  
      typedef UBspline_3d_s       SingleSplineType;  
      typedef BCtype_s            BCType;
      typedef float real_type;
      typedef float value_type;
    };

  template<>
    struct bspline_traits<double,3>
    {
      typedef multi_UBspline_3d_d SplineType;  
      typedef UBspline_3d_d       SingleSplineType;  
      typedef BCtype_d            BCType;
      typedef double real_type;
      typedef double value_type;
    };

  /** specialization for 1D float */
  template<>
    struct bspline_traits<float,1>
    {
      typedef multi_UBspline_1d_s SplineType;
      typedef UBspline_1d_s       SingleSplineType;
      typedef BCtype_s            BCType;
      typedef float               DataType;
      typedef float               real_type;
      typedef float               value_type;
    };

  /** specialization for 1D double */
  template<>
    struct bspline_traits<double,1>
    {
      typedef multi_UBspline_1d_d SplineType;
      typedef UBspline_1d_d       SingleSplineType;
      typedef BCtype_d            BCType;
      typedef double              DataType;
      typedef double              real_type;
      typedef double              value_type;
    };


  /** helper class to determine the value_type of einspline objects
   */
  template<typename ST> struct bspline_type {};

  template<> struct bspline_type<multi_UBspline_3d_s> 
  {
    typedef float value_type;
  };

  template<> struct bspline_type<multi_UBspline_3d_d> 
  {
    typedef double value_type;
  };

  template<> struct bspline_type<UBspline_3d_s> 
  {
    typedef float value_type;
  };

  template<> struct bspline_type<UBspline_3d_d> 
  {
    typedef double value_type;
  };
}
#endif
