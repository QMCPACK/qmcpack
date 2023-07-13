//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from ParticleSet.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <type_traits>
#include "ContextForSteps.h"

namespace qmcplusplus
{
ContextForSteps::ContextForSteps(RandomBase<FullPrecRealType>& random_gen) : random_gen_(random_gen) {}

RandomBase<ContextForSteps::FullPrecRealType>& ContextForSteps::get_random_gen() { return random_gen_; }

} // namespace qmcplusplus
