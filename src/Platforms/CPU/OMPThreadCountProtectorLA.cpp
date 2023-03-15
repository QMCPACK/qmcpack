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

#include "OMPThreadCountProtectorLA.hpp"

namespace qmcplusplus
{
  OMPThreadCountProtectorLA::OMPThreadCountProtectorLA()
  {
#if defined(HAVE_OPENBLAS)
    handle_ = std::make_unique<Protector>();
#endif
  }

  OMPThreadCountProtectorLA::~OMPThreadCountProtectorLA() = default;
} // namespace qmcplusplus
