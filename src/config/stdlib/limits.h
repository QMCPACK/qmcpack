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


#ifndef QMCPLUSPLUS_LIMITS_H
#define QMCPLUSPLUS_LIMITS_H
#include <limits>
#include <type_traits>

#ifndef HAVE_ISZERO
template<typename T,
  typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
inline bool iszero(T a)
{
  return (std::abs(a)<std::numeric_limits<T>::epsilon());
}
#endif

#endif

