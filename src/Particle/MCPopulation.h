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

#include <iterator>
#include <vector>
#include <memory>

#include "Configuration.h"
#include "type_traits/TypeRequire.hpp"
#include "Particle/ParticleSet.h"
#include "ParticleBase/ParticleAttrib.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/Walker.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
class MCPopulation
{
public:
  using MCPWalker  = Walker<QMCTraits, PtclOnLatticeTraits>;
  using RealType   = QMCTraits::RealType;
  using Properties = MCPWalker::PropertyContainer_t;
  using IndexType  = QMCTraits::IndexType;

private:
  int num_ranks_                = 0;
  IndexType num_global_walkers_ = 0;
  IndexType num_local_walkers_  = 0;
  IndexType num_particles_      = 0;
  IndexType num_groups_         = 0;
  IndexType max_samples_        = 0;
  IndexType target_population_  = 0;
  IndexType target_samples_     = 0;
  //Properties properties_;
  ParticleSet ions_;
  std::vector<IndexType> walker_offsets_;
  // By making this a linked list and creating the crowds at the same time we could get first touch.
  std::vector<std::unique_ptr<MCPWalker>> walkers_;
  std::vector<std::pair<int, int>> particle_group_indexes_;
  SpeciesSet species_set_;
  std::vector<RealType> ptclgrp_mass_;
  ///1/Mass per species
  std::vector<RealType> ptclgrp_inv_mass_;
  ///1/Mass per particle
  std::vector<RealType> ptcl_inv_mass_;

  std::shared_ptr<TrialWaveFunction> trial_wf_;
  std::shared_ptr<ParticleSet> elec_particle_set_;
  // At the moment these are "clones" but I think this design pattern smells.
  std::vector<std::unique_ptr<ParticleSet>> walker_elec_particle_sets_;
  std::vector<std::unique_ptr<TrialWaveFunction>> walker_trial_wavefunctions_;
  std::vector<std::unique_ptr<QMCHamiltonian>> walker_hamiltonians_;

public:
  MCPopulation(){};
  MCPopulation(int num_ranks) : num_ranks_(num_ranks) {}
  MCPopulation(MCWalkerConfiguration& mcwc);
  MCPopulation(int num_ranks, int num_particles) : num_ranks_(num_ranks), num_particles_(num_particles) {}
  void createWalkers();
  void createWalkers(IndexType num_walkers, const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& pos);
  void createWalkers(int num_crowds_,
                     int num_walkers_per_crowd_,
                     IndexType num_walkers,
                     const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& pos);

  /** puts walkers and their "cloned" things into groups in a somewhat general way
   *
   *  Should compile only if ITER is a proper input ITERATOR
   *  and implements addWalker
   *  
   */
  template<typename ITER, typename = RequireInputIterator<ITER>>
  void distributeWalkers(ITER it_group, ITER group_end, int walkers_per_group)
  {
    auto it_walkers             = walkers_.begin();
    auto it_walker_elecs        = walker_elec_particle_sets_.begin();
    auto it_walker_twfs         = walker_trial_wavefunctions_.begin();
    auto it_walker_hamiltonians = walker_hamiltonians_.begin();

    while (it_group != group_end)
    {
      for (int i = 0; i < walkers_per_group; ++i)
      {
        // possible that walkers_all < walkers_per_group * group_size
        if (it_walkers == walkers_.end())
          break;
        (**it_group).addWalker(**it_walkers, **it_walker_elecs, **it_walker_twfs, **it_walker_hamiltonians);
        ++it_walkers;
        ++it_walker_elecs;
        ++it_walker_twfs;
        ++it_walker_hamiltonians;
      }
      ++it_group;
    }
  }
  /**@ingroup Accessors
   * @{
   */
  IndexType get_active_walkers() { return walkers_.size(); }
  int get_num_ranks() const { return num_ranks_; }
  IndexType get_num_global_walkers() const { return num_global_walkers_; }
  IndexType get_num_local_walkers() const { return num_local_walkers_; }
  IndexType get_num_particles() const { return num_particles_; }
  IndexType get_max_samples() const { return max_samples_; }
  IndexType get_target_population() const { return target_population_; }
  IndexType get_target_samples() const { return target_samples_; }
  //const Properties& get_properties() const { return properties_; }
  const SpeciesSet& get_species_set() const { return species_set_; }
  const ParticleSet& get_ions() const { return ions_; }
  const std::vector<int>& get_walker_offsets() const { return walker_offsets_; }

  void set_num_global_walkers(IndexType num_global_walkers) { num_global_walkers_ = num_global_walkers; }
  void set_num_local_walkers(IndexType num_local_walkers) { num_local_walkers_ = num_local_walkers; }

  void set_target(IndexType pop) { target_population_ = pop; }
  void set_target_samples(IndexType samples) { target_samples_ = samples; }

  std::vector<std::unique_ptr<MCPWalker>>& get_walkers() { return walkers_; }
  const std::vector<std::pair<int, int>>& get_particle_group_indexes() const { return particle_group_indexes_; }
  const std::vector<RealType>& get_ptclgrp_mass() const { return ptclgrp_mass_; }
  const std::vector<RealType>& get_ptclgrp_inv_mass() const { return ptclgrp_inv_mass_; }
  const std::vector<RealType>& get_ptcl_inv_mass() const { return ptcl_inv_mass_; }

  /** }@ */
};


} // namespace qmcplusplus

#endif
