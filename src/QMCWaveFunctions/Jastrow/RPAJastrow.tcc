//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 Jeongnim Kim and QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_RPA_JASTROW_TCC
#define QMCPLUSPLUS_RPA_JASTROW_TCC

#include "RPAJastrow.h"


namespace qmcplusplus
{

template<class T>
T RPAJastrow::do_ratioT(ParticleSet& P, int iat)
{
  T r(1.0);
  for (int i = 0; i < Psi.size(); i++)
    r *= Psi[i]->ratio<T>(P, iat);
  return r;
}


} // namespace qmcplusplus

#endif