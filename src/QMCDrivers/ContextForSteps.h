//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from ParticleSet.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MOVECONTEXT_H
#define QMCPLUSPLUS_MOVECONTEXT_H

#include <map>
#include <vector>
#include "OhmmsSoA/VectorSoaContainer.h"
#include "Configuration.h"
#include "Particle/Walker.h"
#include "QMCDrivers/Crowd.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus
{
/** Thread local context for moving walkers
 *
 *  created once per driver per crowd
 *  It's two significant responsibilities are holding the thread local RandomGen_t
 *  And the particle group indexes.
 *
 *  
 */
class ContextForSteps
{
public:
  ContextForSteps(RandomGenerator& random_gen);

  RandomGenerator& get_random_gen() { return random_gen_; }

protected:
  RandomGenerator& random_gen_;
};

} // namespace qmcplusplus
#endif
