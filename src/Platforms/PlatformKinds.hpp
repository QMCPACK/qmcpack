//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_PLATFORM_KINDS_H
#define QMCPLUSPLUS_PLATFORM_KINDS_H

namespace qmcplusplus
{

enum class PlatformKind
{
  CPU,
  OMPTARGET,
  CUDA,
  SYCL
};

} // namespace qmcplusplus
#endif
