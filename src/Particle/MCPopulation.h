//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from VMC.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MCPOPULATION_H
#define QMCPLUSPLUS_MCPOPULATION_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/Walker.h"

namespace qmcplusplus
{

class MCPopulation
{
public:
  using PopulationWalker = Walker<QMCTraits, PtclOnLatticeTraits>;
  using Properties = PopulationWalker::PropertyContainer_t;
  using IndexType = QMCTraits::IndexType;
private:
  IndexType global_num_walkers_;
  IndexType local_num_walkers_;
  IndexType max_samples_;
  IndexType target_population_;
  IndexType target_samples_;
  Properties properties_;
  ParticleSet ions_;
  std::vector<IndexType> walker_offsets_;
public:
  MCPopulation() {};
  MCPopulation(MCWalkerConfiguration& mcwc);
  const ParticleSet& get_ions() const {return ions_;}
  const IndexType get_target() const { return target_population_; }
  
  const IndexType get_max_samples() const { return max_samples_; }
  const Properties& get_properties() const {return properties_; }
  const IndexType get_num_global_walkers() const {return global_num_walkers_; }
  const std::vector<int>& get_walker_offsets() const {return walker_offsets_; }
  void set_target(IndexType pop) { target_population_ = pop; }
  void set_target_samples(IndexType samples) { target_samples_ = samples; }
};


}

#endif
