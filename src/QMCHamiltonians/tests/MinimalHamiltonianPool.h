//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MINIMALHAMILTONIANPOOL_H
#define QMCPLUSPLUS_MINIMALHAMILTONIANPOOL_H

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "Particle/ParticleSetPool.h"

namespace qmcplusplus
{
class MinimalHamiltonianPool
{
public:
  static HamiltonianPool make_hamWithEE(Communicate* comm,
                                        ParticleSetPool& particle_pool,
                                        WaveFunctionPool& wavefunction_pool);
  static HamiltonianPool makeHamWithEEEI(Communicate* comm,
                                         ParticleSetPool& particle_pool,
                                         WaveFunctionPool& wavefunction_pool);
};

} // namespace qmcplusplus
#endif
