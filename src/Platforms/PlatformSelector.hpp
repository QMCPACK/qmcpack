//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_PLATFORM_SELECTOR_H
#define QMCPLUSPLUS_PLATFORM_SELECTOR_H

#include <vector>
#include <string>

namespace qmcplusplus
{

enum class PlatformKind
{
  CPU,
  OMPTARGET,
  CUDA,
  SYCL
};

enum class SelectorKind
{
  CPU_OMPTARGET,
  CPU_OMPTARGET_CUDA,
  CPU_OMPTARGET_SYCL,
};

template<SelectorKind KIND>
class PlatformSelector
{
public:
  static const std::vector<std::string> candidate_values;
  static PlatformKind selectPlatform(std::string_view value);
};

using CPUOMPTargetSelector = PlatformSelector<SelectorKind::CPU_OMPTARGET>;

template<>
const std::vector<std::string> PlatformSelector<SelectorKind::CPU_OMPTARGET>::candidate_values;
template<>
const std::vector<std::string> PlatformSelector<SelectorKind::CPU_OMPTARGET_CUDA>::candidate_values;
template<>
const std::vector<std::string> PlatformSelector<SelectorKind::CPU_OMPTARGET_SYCL>::candidate_values;

#if defined(ENABLE_CUDA)
using CPUOMPTargetVendorSelector = PlatformSelector<SelectorKind::CPU_OMPTARGET_CUDA>;
#elif defined(ENABLE_SYCL)
using CPUOMPTargetVendorSelector = PlatformSelector<SelectorKind::CPU_OMPTARGET_SYCL>;
#else
using CPUOMPTargetVendorSelector = PlatformSelector<SelectorKind::CPU_OMPTARGET>;
#endif

} // namespace qmcplusplus
#endif
