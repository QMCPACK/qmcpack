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
#ifndef QMCPLUSPLUS_RANDOM_FORWARD_H
#define QMCPLUSPLUS_RANDOM_FORWARD_H

inline uint32_t MakeSeed(int i, int n)
{
  const uint32_t u=1<<10;
  return static_cast<uint32_t>(std::time(nullptr))%u+(i+1)*n+i;
}

#include <Utilities/RandomGenerator.h>
#endif

