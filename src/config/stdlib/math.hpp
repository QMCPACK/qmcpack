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
#if defined(HAVE_MASS)
#include <mass.h>
#endif

#if __APPLE__

inline void sincos(double a, double* s, double* c)
{
  __sincos(a,s,c);
}

inline void sincos(float a, float* s, float* c)
{
  __sincosf(a,s,c);
}

#else // not __APPLE__

#if defined(HAVE_SINCOS)

inline void sincos(float a, float* s, float* c)
{
#if defined(HAVE_MASS)
  // there is no sincosf in libmass
  // libmass sincos is faster than libm sincosf
  double ds,dc;
  sincos((double)a,&ds,&dc);
  *s=ds; *c=dc;
#else
  sincosf(a,s,c);
#endif
}

#else

template<typename T>
inline void sincos(T a, T* restrict s, T*  restrict c)
{
  *s=std::sin(a);
  *c=std::cos(a);
}

#endif // HAVE_SINCOS

#endif // __APPLE__

namespace qmcplusplus
{
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

}

#endif
