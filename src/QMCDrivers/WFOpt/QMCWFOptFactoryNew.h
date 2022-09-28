//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_WFOPTFACTORYNEW_H
#define QMCPLUSPLUS_WFOPTFACTORYNEW_H

#include "QMCDrivers/QMCDriverInterface.h"

class Communicate;

namespace qmcplusplus
{
class MCPopulation;
class WaveFunctionPool;
class QMCHamiltonian;
class TrialWaveFunction;
class SampleStack;
class QMCFixedSampleLinearOptimizeBatched;
class ProjectData;

std::unique_ptr<QMCFixedSampleLinearOptimizeBatched> QMCWFOptLinearFactoryNew(
    xmlNodePtr cur,
    const ProjectData& project_data,
    const std::optional<EstimatorManagerInput>& global_emi,
    WalkerConfigurations& wc,
    MCPopulation&& pop,
    SampleStack& samples,
    Communicate* comm);
} // namespace qmcplusplus

#endif
