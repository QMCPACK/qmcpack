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

/// 1/2 - 1/2 tanh(alpha * (x - 1/2))
template<typename T>
T smooth_func_tanh(T x, T& dx, T& d2x)
{
  if(x<0)
  {
    dx = d2x = T(0);
    return T(1);
  }
  else if(x>=1)
  {
    dx = d2x = T(0);
    return T(0);
  }
  else
  {
    const T cone(1), chalf(0.5), alpha(2);
    const T tanh_x = std::tanh((x - chalf) * alpha);
    const T dtanhx_dx = cone - tanh_x * tanh_x;

    dx = - chalf * alpha * dtanhx_dx;
    d2x = alpha * alpha * tanh_x * dtanhx_dx;
    return chalf * (cone - tanh_x);
  }
}

} // namespace qmcplusplus

#if defined(DEBUG_MAIN)
#include <iostream>

int main()
{
  using RealType = double;
  const RealType num_grid(1000);
  std::cout << "# x      dx     d2x" << std::endl;
  for(int i=-2; i<num_grid+2; i++)
  {
    RealType x = RealType(i)/num_grid;
    RealType f, df, d2f;
    f = qmcplusplus::smooth_func_tanh(x, df, d2f);
    std::cout << x << " " << f << " " << df << " " << d2f << std::endl;
  }
  return 0;
}
#endif
#endif
