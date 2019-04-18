//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SMOOTH_FUNCTIONS_H
#define QMCPLUSPLUS_SMOOTH_FUNCTIONS_H

#include <cmath>

namespace qmcplusplus
{
struct SmoothFunctions
{
  enum
  {
    YE2018 = 0,
    COSCOS,
    LINEAR
  };

  /// 1/2 - 1/2 tanh(alpha * (x - 1/2))
  template<typename T>
  static T func_tanh(T x, T& dx, T& d2x)
  {
    if (x < 0)
    {
      dx = d2x = T(0);
      return T(1);
    }
    else if (x >= 1)
    {
      dx = d2x = T(0);
      return T(0);
    }
    else
    {
      const T cone(1), chalf(0.5), alpha(2);
      const T tanh_x    = std::tanh((x - chalf) * alpha);
      const T dtanhx_dx = cone - tanh_x * tanh_x;

      dx  = -chalf * alpha * dtanhx_dx;
      d2x = alpha * alpha * tanh_x * dtanhx_dx;
      return chalf * (cone - tanh_x);
    }
  }

  /// (1+cos(PI*(1-cos(PI*x))/2))/2
  template<typename T>
  static T func_coscos(T x, T& dx, T& d2x)
  {
    if (x < 0)
    {
      dx = d2x = T(0);
      return T(1);
    }
    else if (x >= 1)
    {
      dx = d2x = T(0);
      return T(0);
    }
    else
    {
      const T chalf(0.5), cone(1), pihalf(M_PI * chalf), pipihalf(M_PI * M_PI * chalf);
      T s, c, scos, ccos;
      sincos(T(M_PI) * x, &s, &c);
      sincos(pihalf * (cone - c), &scos, &ccos);

      dx  = -chalf * pipihalf * scos * s;
      d2x = -pihalf * pipihalf * (ccos * pihalf * s * s + scos * c);
      return chalf * (cone + ccos);
    }
  }

  /// 1-x
  template<typename T>
  static T func_linear(T x, T& dx, T& d2x)
  {
    if (x < 0)
    {
      dx = d2x = T(0);
      return T(1);
    }
    else if (x >= 1)
    {
      dx = d2x = T(0);
      return T(0);
    }
    else
    {
      dx  = T(-1);
      d2x = T(0);
      return T(1) - x;
    }
  }
};

} // namespace qmcplusplus

#if defined(DEBUG_MAIN)
#include <iostream>

int main()
{
  using RealType = double;
  const RealType num_grid(1000);
  std::cout << "# x      dx     d2x" << std::endl;
  for (int i = -2; i < num_grid + 2; i++)
  {
    RealType x = RealType(i) / num_grid;
    RealType f, df, d2f;
    f = qmcplusplus::SmoothFunctions::func_coscos(x, df, d2f);
    std::cout << x << " " << f << " " << df << " " << d2f << std::endl;
  }
  return 0;
}
#endif
#endif
