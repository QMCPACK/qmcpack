//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_STDLIB_PORT_H
#define QMCPLUSPLUS_STDLIB_PORT_H
#include <config.h>
#include <cmath>
#include <type_traits>
#include "config/stdlib/Constants.h"
#if defined(HAVE_AMD_LIBM)
#include <amdlibm.h>
#endif
#if defined(HAVE_MASS)
#include <mass.h>
#endif

namespace qmcplusplus
{

/// sincos function wrapper
#if defined(__APPLE__)

inline void sincos(double a, double* restrict s, double* restrict c)
{
  ::__sincos(a,s,c);
}

inline void sincos(float a, float* restrict s, float* restrict c)
{
  ::__sincosf(a,s,c);
}

#elif defined(HAVE_AMD_LIBM)

inline void sincos(double a, double* restrict s, double* restrict c)
{
  ::amd_sincos(a,s,c);
}

inline void sincos(float a, float* restrict s, float* restrict c)
{
  ::amd_sincosf(a,s,c);
}

#elif defined(HAVE_SINCOS)

inline void sincos(double a, double* restrict s, double* restrict c)
{
  ::sincos(a,s,c);
}

inline void sincos(float a, float* restrict s, float* restrict c)
{
#if defined(HAVE_MASS)
  // there is no sincosf in libmass
  // libmass sincos is faster than libm sincosf
  double ds,dc;
  ::sincos((double)a,&ds,&dc);
  *s=ds; *c=dc;
#else
  ::sincosf(a,s,c);
#endif
}

#else // fallback

template<typename T>
inline void sincos(T a, T* restrict s, T*  restrict c)
{
  *s=std::sin(a);
  *c=std::cos(a);
}

#endif

/** return i^n
 *
 * std::pow(int,int) is not standard
 */
inline int pow(int i, int n)
{
  return static_cast<int>(std::pow(static_cast<double>(i),n));
}

template<typename T,
  typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
inline bool iszero(T a)
{
  return std::fpclassify(a) == FP_ZERO;
}

/** return true if the value is NaN.
 * std::isnan can be affected by compiler the -ffast-math option and return true constantly.
 * The customized qmcplusplus::isnan should be always effective.
 * This requires its definition compiled without -ffast-math.
 */
bool isnan(float);
bool isnan(double);

}

#endif
