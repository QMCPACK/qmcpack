//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OMPTARGET_MATH_H
#define OMPTARGET_MATH_H

#include <cmath>

namespace omptarget
{
inline void sincos(double a, double* restrict s, double* restrict c)
{
  ::sincos(a,s,c);
}

inline void sincos(float a, float* restrict s, float* restrict c)
{
  ::sincosf(a,s,c);
}
}
#endif
