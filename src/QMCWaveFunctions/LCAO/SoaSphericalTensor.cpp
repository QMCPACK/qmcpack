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


#include "SoaSphericalTensor.h"

namespace qmcplusplus
{
template struct SoaSphericalTensor<float>;
template struct SoaSphericalTensor<double>;
} // namespace qmcplusplus
