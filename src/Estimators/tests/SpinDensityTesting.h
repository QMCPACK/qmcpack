//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPINDENSITYTESTING_H
#define QMCPLUSPLUS_SPINDENSITYTESTING_H

#include "ParticleSet.h"

namespace qmcplusplus
{
namespace testing
{

using POLT    = PtclOnLatticeTraits;
using Lattice = POLT::ParticleLayout_t;

Lattice makeTestLattice();

}
}
#endif
