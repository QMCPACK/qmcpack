//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MINIMALWAVEFUNCTIONPOOL_H
#define QMCPLUSPLUS_MINIMALWAVEFUNCTIONPOOL_H

#include "Message/Communicate.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"

namespace qmcplusplus
{

class MinimalWaveFunctionPool
{
public:
  static WaveFunctionPool make_diamondC_1x1x1(const RuntimeOptions& runtime_options,
                                              Communicate* comm,
                                              ParticleSetPool& particle_pool);
  static WaveFunctionPool make_O2_spinor(const RuntimeOptions& runtime_options,
                                         Communicate* comm,
                                         ParticleSetPool& particle_pool);
  static WaveFunctionPool make_O2_spinor_J12(const RuntimeOptions& runtime_options,
                                             Communicate* comm,
                                             ParticleSetPool& particle_pool);
};

} // namespace qmcplusplus
#endif
