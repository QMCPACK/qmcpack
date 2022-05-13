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

template<>
const std::vector<std::string> PlatformSelector<SelectorKind::CPU_OMPTARGET>::candidate_values{"", "yes", "no",
                                                                                               "omptarget", "cpu"};

template<>
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

template<>
const std::vector<std::string> PlatformSelector<SelectorKind::CPU_OMPTARGET_CUDA>::candidate_values{"",          "yes",
                                                                                                    "no",        "cuda",
                                                                                                    "omptarget", "cpu"};

template<>
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

template<>
const std::vector<std::string> PlatformSelector<SelectorKind::CPU_OMPTARGET_SYCL>::candidate_values{"",          "yes",
                                                                                                    "no",        "sycl",
                                                                                                    "omptarget", "cpu"};

template<>
PlatformKind PlatformSelector<SelectorKind::CPU_OMPTARGET_SYCL>::selectPlatform(std::string_view value)
{
#if defined(ENABLE_SYCL)
  if (value.empty() || value == "yes" || value == "sycl")
    return PlatformKind::SYCL;
  else if (value == "omptarget")
#else
  if (value == "sycl")
    throw std::runtime_error("Cannot access SYCL code path. Executables are built with ENABLE_SYCL=OFF.");
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
