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

/** @file MCSample.cpp
 * @brief Stores particle configurations for later use in DMC and wavefunction optimization
 *
 * Stores less temporary data than the buffer object.
 */

#include "Particle/MCSample.h"

namespace qmcplusplus
{
using Walker_t = ParticleSet::Walker_t;

MCSample convertWalkerToSample(Walker_t& walker)
{
  int n = walker.R.size();
  MCSample sample(n);
  sample.put(walker);
  return sample;
}


} // namespace qmcplusplus
