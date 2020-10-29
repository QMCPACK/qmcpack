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

namespace qmcplusplus
{

enum class smoothing_functions
{
  LEKS2018 = 0,
  COSCOS,
  LINEAR
};

template<typename T>
T smoothing(smoothing_functions func_id, T x, T& dx, T& d2x);

extern template float smoothing(smoothing_functions func_id, float x, float& dx, float& d2x);
extern template double smoothing(smoothing_functions func_id, double x, double& dx, double& d2x);
}

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
    f = qmcplusplus::SmoothFunctions::func(COSCOS, x, df, d2f);
    std::cout << x << " " << f << " " << df << " " << d2f << std::endl;
  }
  return 0;
}
#endif
#endif
