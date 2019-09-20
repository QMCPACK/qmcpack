//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from DMC.h
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DMCBATCHED_H
#define QMCPLUSPLUS_DMCBATCHED_H

#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/DMC/DMCDriverInput.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/ContextForSteps.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a DMC using particle-by-particle move. Threaded execution.
 */
class DMCBatched : public QMCDriverNew
{
public:
  using FullPrecRealType  = QMCTraits::FullPrecRealType;
  using PosType           = QMCTraits::PosType;
  using ParticlePositions = PtclOnLatticeTraits::ParticlePos_t;

  /** To avoid 10's of arguments to runDMCStep
   *
   *  There should be a division between const input to runVMCStep
   *  And step to step state
   */
  struct StateForThread
  {
    const QMCDriverInput& qmcdrv_input;
    const DMCDriverInput& dmcdrv_input;
    const MCPopulation& population;
    IndexType recalculate_properties_period;
    IndexType step;
    int block;
    bool recomputing_blocks;
    const DriftModifierBase& drift_modifier;
    StateForThread(QMCDriverInput& qmci, DMCDriverInput& vmci, DriftModifierBase& drift_mod, MCPopulation& pop)
        : qmcdrv_input(qmci), dmcdrv_input(vmci), drift_modifier(drift_mod), population(pop)
    {}
  };

  /// Constructor.
  DMCBatched(QMCDriverInput&& qmcdriver_input,
             DMCDriverInput&& input,
             MCPopulation&& pop,
             TrialWaveFunction& psi,
             QMCHamiltonian& h,
             WaveFunctionPool& ppool,
             Communicate* comm);

  bool run();

private:
  DMCDriverInput dmcdriver_input_;
  ///Interval between branching
  IndexType branch_interval_;
  void resetUpdateEngines();
  /// Copy Constructor (disabled)
  DMCBatched(const DMCBatched&) = delete;
  /// Copy operator (disabled).
  DMCBatched& operator=(const DMCBatched&) = delete;
};

} // namespace qmcplusplus

#endif
