//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**@file ExampleHeBuilder.h
 *@brief Example builder for simple He wavefunction.
 */
#ifndef QMCPLUSPLUS_EXAMPLEHEBUILDER_H
#define QMCPLUSPLUS_EXAMPLEHEBUILDER_H

#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

/**@defgroup WFSBuilder Orbital builder group
 * @brief Builder classes to add WaveFunctionComponent to a TrialWaveFunction
 */
namespace qmcplusplus
{
class ExampleHeBuilder : public WaveFunctionComponentBuilder
{
public:
  ExampleHeBuilder(Communicate* comm, ParticleSet& p, const PtclPoolType& psets);

  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;

private:
  const PtclPoolType& ptclPool;
  ParticleSet& els;
};


} // namespace qmcplusplus
#endif
