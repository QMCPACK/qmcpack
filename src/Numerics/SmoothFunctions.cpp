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

#include <stdexcept>

#include "CPU/math.hpp"
#include "SmoothFunctions.hpp"

namespace qmcplusplus
{

template<typename T>
T smoothing(smoothing_functions func_id, T x, T& dx, T& d2x)
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
  else if (func_id == smoothing_functions::LEKS2018)
  {
    /// 1/2 - 1/2 tanh(alpha * (x - 1/2))
    const T cone(1), chalf(0.5), alpha(2);
    const T tanh_x    = std::tanh((x - chalf) * alpha);
    const T dtanhx_dx = cone - tanh_x * tanh_x;

    dx  = -chalf * alpha * dtanhx_dx;
    d2x = alpha * alpha * tanh_x * dtanhx_dx;
    return chalf * (cone - tanh_x);
  }
  else if (func_id == smoothing_functions::COSCOS)
  {
    /// (1+cos(PI*(1-cos(PI*x))/2))/2
    const T chalf(0.5), cone(1), pihalf(M_PI * chalf), pipihalf(M_PI * M_PI * chalf);
    T s, c, scos, ccos;
    qmcplusplus::sincos(T(M_PI) * x, &s, &c);
    qmcplusplus::sincos(pihalf * (cone - c), &scos, &ccos);

    dx  = -chalf * pipihalf * scos * s;
    d2x = -pihalf * pipihalf * (ccos * pihalf * s * s + scos * c);
    return chalf * (cone + ccos);
  }
  else if (func_id == smoothing_functions::LINEAR)
  {
    /// 1-x
    dx  = T(-1);
    d2x = T(0);
    return T(1) - x;
  }
  else
    throw std::runtime_error("Unknown smooth function!");
}

template float smoothing(smoothing_functions func_id, float x, float& dx, float& d2x);
template double smoothing(smoothing_functions func_id, double x, double& dx, double& d2x);

} // namespace qmcplusplus


