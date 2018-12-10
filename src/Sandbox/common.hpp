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
#ifndef QMCPLUSPLUS_MINIAPPS_COMMON_H
#define QMCPLUSPLUS_MINIAPPS_COMMON_H
/** @file common.hpp
 * @brief Testing random number generators
 */
#include <Utilities/Timer.h>
#include <random/random.hpp>
#include <vector>
#include <getopt.h>

namespace qmcplusplus
{
  struct Square
  {
    template<typename T>
      inline T operator()(const T &lhs, const T& rhs) const 
      {
        return lhs+rhs*rhs;
      }
  };

  template<typename T>
    struct Stat
    {
      T sum;
      T sumsq;

      Stat():sum(T()),sumsq(T()) {}

      inline typename std::pair<T,T> apply(const T* d, int n)
      {
        sum=std::accumulate(d,d+n,T());
        sumsq=std::accumulate(d,d+n,T(),Square());
        T avg=sum/n;
        return make_pair(avg,sumsq/n-avg*avg);
      }
    };
}
#endif
