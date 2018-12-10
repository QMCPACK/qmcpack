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
#include <einspline/multi_bspline.h>
namespace qmcplusplus
{
  /** determine if EngT (e.g., einspline engine) handles real data or complex data
   *
   * Default is true and complex cases are specialized
   */
  template<typename EngT>
  struct is_real_bspline
  {
      static const bool value=true;
  };

  ///specialization for multi_UBspline_3d_z
  template<>
  struct is_real_bspline<multi_UBspline_3d_z>
  {
      static const bool value=false;
  };

  ///specialization for multi_UBspline_3d_z
  template<>
  struct is_real_bspline<multi_UBspline_3d_c>
  {
      static const bool value=false;
  };
  //
  /** struct to check if two types are the same
   */
  template<typename T1, typename T2>
    struct type_check
    {
      static const bool value=false;
    };

  template<typename T1>
    struct type_check<T1,T1>
    {
      static const bool value=true;
    };

  template<typename T1, typename T2> struct is_same_precision
  {
    static const bool value=false;
  };
  template<> struct is_same_precision<double,double>
  {
    static const bool value=true;
  };
  template<> struct is_same_precision<float,float>
  {
    static const bool value=true;
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
  struct bspline_engine_traits {};

#if OHMMS_DIM == 3
  /** specialization with multi_UBspline_3d_d
   */
  template<>
  struct bspline_engine_traits<multi_UBspline_3d_d>
  {
    enum {DIM=3};
    typedef multi_UBspline_3d_d SplineType;  
    typedef UBspline_3d_d       SingleSplineType;  
    typedef BCtype_d            BCType;
    typedef double real_type;
    typedef double value_type;
  };

  ///specialization with multi_UBspline_3d_z
  template<>
  struct bspline_engine_traits<multi_UBspline_3d_z>
  {
    enum {DIM=3};
    typedef multi_UBspline_3d_z SplineType;  
    typedef UBspline_3d_z       SingleSplineType;  
    typedef BCtype_z            BCType;
    typedef std::complex<double> value_type;
    typedef UBspline_3d_z single_spline_type;
  };

  /** specialization with multi_UBspline_3d_s
   */
  template<>
  struct bspline_engine_traits<multi_UBspline_3d_s>
  {
    enum {DIM=3};
    typedef multi_UBspline_3d_s SplineType;  
    typedef UBspline_3d_s       SingleSplineType;  
    typedef BCtype_s            BCType;
    typedef float real_type;
    typedef float value_type;
  };

  ///specialization with multi_UBspline_3d_c
  template<>
  struct bspline_engine_traits<multi_UBspline_3d_c>
  {
    enum {DIM=3};
    typedef multi_UBspline_3d_c SplineType;  
    typedef UBspline_3d_c       SingleSplineType;  
    typedef BCtype_c            BCType;
    typedef float real_type;
    typedef std::complex<float> value_type;
  };

#else
  /** specialization with multi_UBspline_3d_d
   */
  template<>
  struct bspline_engine_traits<multi_UBspline_2d_d>
  {
    enum {DIM=2};
    typedef double real_type;
    typedef double value_type;
  };

  ///specialization with multi_UBspline_3d_z
  template<>
  struct bspline_engine_traits<multi_UBspline_2d_z>
  {
    enum {DIM=2};
    typedef double real_type;
    typedef std::complex<double> value_type;
  };

  /** specialization with multi_UBspline_3d_d
   */
  template<>
  struct bspline_engine_traits<multi_UBspline_2d_s>
  {
    enum {DIM=2};
    typedef float real_type;
    typedef float value_type;
  };

  ///specialization with multi_UBspline_3d_z
  template<>
  struct bspline_engine_traits<multi_UBspline_2d_c>
  {
    enum {DIM=2};
    typedef float real_type;
    typedef std::complex<float> value_type;
  };
#endif

  /** specialization with multi_UBspline_3d_d
   */
  template<>
  struct bspline_engine_traits<multi_UBspline_1d_d>
  {
    enum {DIM=1};
    typedef multi_UBspline_1d_d SplineType;
    typedef UBspline_1d_d       SingleSplineType;
    typedef BCtype_d            BCType;
    typedef double real_type;
    typedef double value_type;
  };

  /** specialization with multi_UBspline_3d_s
   */
  template<>
  struct bspline_engine_traits<multi_UBspline_1d_s>
  {
    enum {DIM=1};
    typedef multi_UBspline_1d_s SplineType;
    typedef UBspline_1d_s       SingleSplineType;
    typedef BCtype_s            BCType;
    typedef float real_type;
    typedef float value_type;
  };

}
#endif
