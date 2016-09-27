//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.   
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc. 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LIMITS_H
#define QMCPLUSPLUS_LIMITS_H
#include <limits>
template<typename T>
inline bool iszero(T a)
{
  return (std::abs(a)<std::numeric_limits<T>::epsilon());
}
#endif

