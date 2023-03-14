//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#include "LAThreadProtector.hpp"

namespace qmcplusplus
{
  LAThreadProtector::LAThreadProtector()
  {
#if defined(HAVE_OPENBLAS)
    handle_ = std::make_unique<Concurrency::ThreadCountProtector<>>();
#endif
  }

  LAThreadProtector::~LAThreadProtector() = default;
} // namespace qmcplusplus
