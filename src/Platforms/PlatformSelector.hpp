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
  CUDA
};

enum class SelectorKind
{
  CPU_OMPTARGET,
  CPU_OMPTARGET_CUDA
};

template<SelectorKind KIND>
class PlatformSelector;

template<>
class PlatformSelector<SelectorKind::CPU_OMPTARGET>
{
public:
  static const std::vector<std::string> candidate_values;
  static PlatformKind selectPlatform(std::string_view value);
};

using CPUOMPTargetSelector = PlatformSelector<SelectorKind::CPU_OMPTARGET>;

template<>
class PlatformSelector<SelectorKind::CPU_OMPTARGET_CUDA>
{
public:
  static const std::vector<std::string> candidate_values;
  static PlatformKind selectPlatform(std::string_view value);
};

using CPUOMPTargetCUDASelector = PlatformSelector<SelectorKind::CPU_OMPTARGET_CUDA>;
} // namespace qmcplusplus
#endif
