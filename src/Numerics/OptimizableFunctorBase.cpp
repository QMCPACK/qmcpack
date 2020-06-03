//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "OptimizableFunctorBase.h"

void print(OptimizableFunctorBase& func, std::ostream& os)
{
  typedef OptimizableFunctorBase::real_type real_type;
  int n       = 100;
  real_type d = func.cutoff_radius / 100., r = 0;
  real_type u, du;
  for (int i = 0; i < n; ++i)
  {
    u  = func.f(r);
    du = func.df(r);
    os << std::setw(22) << r << std::setw(22) << u << std::setw(22) << du << std::endl;
    r += d;
  }
}
