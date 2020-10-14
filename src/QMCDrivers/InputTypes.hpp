//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCDriver.cpp
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_INPUTTYPES_H
#define QMCPLUSPLUS_INPUTTYPES_H

#include "Configuration.h"

namespace qmcplusplus
{
namespace input
{
struct PeriodStride
{
  using IndexType  = QMCTraits::IndexType;
  IndexType stride = 0;
  IndexType period = 0;
};
} // namespace input
} // namespace qmcplusplus
#endif
