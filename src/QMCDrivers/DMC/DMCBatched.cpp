//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: DMC.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/DMC/DMCBatched.h"
#include "Concurrency/TasksOneToOne.hpp"
#include "Concurrency/Info.hpp"
#include "Utilities/RunTimeManager.h"
#include "ParticleBase/RandomSeqGenerator.h"


namespace qmcplusplus
{
/** Constructor maintains proper ownership of input parameters
   */
DMCBatched::DMCBatched(QMCDriverInput&& qmcdriver_input,
                       DMCDriverInput&& input,
                       MCPopulation&& pop,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       WaveFunctionPool& ppool,
                       Communicate* comm)
    : QMCDriverNew(std::move(qmcdriver_input), std::move(pop), psi, h, ppool, "DMCBatched::", comm),
      dmcdriver_input_(input)
{
  QMCType = "DMCBatched";
}

} // namespace qmcplusplus
