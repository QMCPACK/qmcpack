//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "HybridRepCenterOrbitalsT.h"

namespace qmcplusplus
{
template class AtomicOrbitalsT<float>;
template class AtomicOrbitalsT<double>;
template class HybridRepCenterOrbitalsT<float, float>;
template class HybridRepCenterOrbitalsT<float, double>;
template class HybridRepCenterOrbitalsT<double, float>;
template class HybridRepCenterOrbitalsT<double, double>;
} // namespace qmcplusplus
