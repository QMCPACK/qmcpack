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
  using Properties       = PopulationWalker::PropertyContainer_t;
  using IndexType        = QMCTraits::IndexType;

private:
  int num_ranks_                = 0;
  IndexType num_global_walkers_ = 0;
  IndexType num_local_walkers_  = 0;
  IndexType max_samples_        = 0;
  IndexType target_population_  = 0;
  IndexType target_samples_     = 0;
  Properties properties_;
  ParticleSet ions_;
  std::vector<IndexType> walker_offsets_;

public:
  MCPopulation(){};
  MCPopulation(int num_ranks) : num_ranks_(num_ranks) {}
  MCPopulation(MCWalkerConfiguration& mcwc);
  int get_num_ranks() const { return num_ranks_; }
  IndexType get_num_global_walkers() const { return num_global_walkers_; }
  IndexType get_num_local_walkers() const { return num_local_walkers_; }
  IndexType get_max_samples() const { return max_samples_; }
  IndexType get_target_population() const { return target_population_; }
  IndexType get_target_samples() const { return target_samples_; }
  const Properties& get_properties() const { return properties_; }
  const ParticleSet& get_ions() const { return ions_; }
  const std::vector<int>& get_walker_offsets() const { return walker_offsets_; }

  void set_num_global_walkers(IndexType num_global_walkers) { num_global_walkers_ = num_global_walkers; }
  void set_num_local_walkers(IndexType num_local_walkers) { num_local_walkers_ = num_local_walkers; }

  void set_target(IndexType pop) { target_population_ = pop; }
  void set_target_samples(IndexType samples) { target_samples_ = samples; }
};


} // namespace qmcplusplus

#endif
