//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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
class QMCOptimizeBatched;
class WaveFunctionPool;
class QMCHamiltonian;
class TrialWaveFunction;
class SampleStack;
class MCWalkerConfiguration;
class HamiltonianPool;
class QMCOptimize;
class QMCFixedSampleLinearOptimizeBatched;
class ProjectData;

QMCOptimizeBatched* QMCWFOptFactoryNew(xmlNodePtr cur,
                                       const ProjectData& project_data,
                                       MCWalkerConfiguration& w,
                                       MCPopulation&& pop,
                                       SampleStack& samples,
                                       Communicate* comm);

QMCFixedSampleLinearOptimizeBatched* QMCWFOptLinearFactoryNew(xmlNodePtr cur,
                                                              const ProjectData& project_data,
                                                              MCWalkerConfiguration& w,
                                                              MCPopulation&& pop,
                                                              SampleStack& samples,
                                                              Communicate* comm);
} // namespace qmcplusplus

#endif
