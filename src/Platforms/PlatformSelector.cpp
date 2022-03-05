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


#include "PlatformSelector.hpp"
#include <config.h>
#include <stdexcept>

namespace qmcplusplus
{

const std::vector<std::string> PlatformSelector<SelectorKind::CPU_OMPTARGET>::candidate_values{"", "yes", "no",
                                                                                               "omptarget", "cpu"};

PlatformKind PlatformSelector<SelectorKind::CPU_OMPTARGET>::selectPlatform(std::string_view value)
{
#if defined(ENABLE_OFFLOAD)
  if (value.empty() || value == "yes" || value == "omptarget")
#else
  if (value == "omptarget")
#endif
    return PlatformKind::OMPTARGET;
  else
    return PlatformKind::CPU;
}

const std::vector<std::string> PlatformSelector<SelectorKind::CPU_OMPTARGET_CUDA>::candidate_values{"",          "yes",
                                                                                                    "no",        "cuda",
                                                                                                    "omptarget", "cpu"};

PlatformKind PlatformSelector<SelectorKind::CPU_OMPTARGET_CUDA>::selectPlatform(std::string_view value)
{
#if defined(ENABLE_CUDA)
  if (value.empty() || value == "yes" || value == "cuda")
    return PlatformKind::CUDA;
  else if (value == "omptarget")
#else
  if (value == "cuda")
    throw std::runtime_error("Cannot access CUDA code path. Executables are built with ENABLE_CUDA=OFF.");
#if defined(ENABLE_OFFLOAD)
  if (value.empty() || value == "yes" || value == "omptarget")
#else
  if (value == "omptarget")
#endif
#endif
    return PlatformKind::OMPTARGET;
  else
    return PlatformKind::CPU;
}

} // namespace qmcplusplus
