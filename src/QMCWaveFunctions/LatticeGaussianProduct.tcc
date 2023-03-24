//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_LATTICE_GAUSSIAN_PRODUCT_TCC
#define QMCPLUSPLUS_LATTICE_GAUSSIAN_PRODUCT_TCC

#include "LatticeGaussianProduct.h"

namespace qmcplusplus
{


template<class T>
T LatticeGaussianProduct::do_ratioT(ParticleSet& P, int iat)
{
  const auto& d_table = P.getDistTableAB(myTableID);
  int icent           = ParticleCenter[iat];
  if (icent == -1)
    return 1.0;
  RealType newdist = d_table.getTempDists()[icent];
  curVal           = ParticleAlpha[iat] * (newdist * newdist);
  return std::exp(static_cast<T>(U[iat] - curVal));
}

} // namespace qmcplusplus

#endif