//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020  QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

/** @file SampleStack.h
 * @brief Stores particle configurations for later use in DMC and wavefunction optimization
 *
 * Stores less temporary data than the buffer object.
 */

#ifndef QMCPLUSPLUS_SAMPLE_STACK_H
#define QMCPLUSPLUS_SAMPLE_STACK_H

#include "Configuration.h"
#include "Particle/SampleStackT.h"

namespace qmcplusplus
{
using SampleStack = SampleStackT<QMCTraits::ValueType>;

} // namespace qmcplusplus
#endif
