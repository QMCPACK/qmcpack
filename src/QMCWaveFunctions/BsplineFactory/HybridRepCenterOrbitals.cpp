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


#include "HybridRepCenterOrbitals.h"

namespace qmcplusplus
{
template class AtomicOrbitals<float>;
template class AtomicOrbitals<double>;
template class HybridRepCenterOrbitals<float>;
template class HybridRepCenterOrbitals<double>;
} // namespace qmcplusplus
