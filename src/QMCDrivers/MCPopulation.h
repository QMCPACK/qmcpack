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

namespace qmcplusplus
{

class QMCHamiltonian;

class MCPopulation
{
public:
  using MCPWalker  = Walker<QMCTraits, PtclOnLatticeTraits>;
  using WFBuffer   = MCPWalker::WFBuffer_t;
  using RealType   = QMCTraits::RealType;
  using Properties = MCPWalker::PropertyContainer_t;
  using IndexType  = QMCTraits::IndexType;

private:
  // Potential thread safety issue
  MCDataType<QMCTraits::FullPrecRealType> ensemble_property_;

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
  std::vector<IndexType> num_local_walkers_per_node_;
  
  // By making this a linked list and creating the crowds at the same time we could get first touch.
  UPtrVector<MCPWalker> walkers_;
  UPtrVector<MCPWalker> dead_walkers_;
  std::vector<std::pair<int, int>> particle_group_indexes_;
  SpeciesSet species_set_;
  std::vector<RealType> ptclgrp_mass_;
  ///1/Mass per species
  std::vector<RealType> ptclgrp_inv_mass_;
  ///1/Mass per particle
  std::vector<RealType> ptcl_inv_mass_;
  size_t size_dataset_;
  // could be
  // std::shared_ptr<TrialWaveFunction> trial_wf_;
  // std::shared_ptr<ParticleSet> elec_particle_set_;
  // std::shared_ptr<QMCHamiltonian> hamiltonian_;

  // This is necessary MCPopulation is constructed in a simple call scope in QMCDriverFactory from the legacy MCWalkerConfiguration
  // MCPopulation should have QMCMain scope eventually and the driver will just have a refrence to it.
  TrialWaveFunction* trial_wf_;
  ParticleSet* elec_particle_set_;
  QMCHamiltonian* hamiltonian_;
  // At the moment these are "clones" but I think this design pattern smells.
  UPtrVector<ParticleSet> walker_elec_particle_sets_;
  UPtrVector<TrialWaveFunction> walker_trial_wavefunctions_;
  UPtrVector<QMCHamiltonian> walker_hamiltonians_;
  // This should perhaps just be aquired from comm but currently MCPopulation
  // is innocent of comm, Every object needing a copy is suboptimal.
  int rank_;

public:
  MCPopulation();
  /** Temporary constructor to deal with MCWalkerConfiguration be the only source of some information
   *  in QMCDriverFactory.
   */
  MCPopulation(int num_ranks,
               MCWalkerConfiguration& mcwc,
               ParticleSet* elecs,
               TrialWaveFunction* trial_wf,
               QMCHamiltonian* hamiltonian_,
               int this_rank = 0);

  MCPopulation(int num_ranks,
               ParticleSet* elecs,
               TrialWaveFunction* trial_wf,
               QMCHamiltonian* hamiltonian,
               int this_rank = 0);

  MCPopulation(MCPopulation&)  = default;
  MCPopulation(MCPopulation&&) = default;
  MCPopulation& operator=(MCPopulation&&) = default;

  /** @ingroup PopulationControl
   *
   *  State Requirement:
   *   * createWalkers must have been called
   *  @{
   */
  MCPWalker& spawnWalker();
  void killWalker(MCPWalker&);
  void killLastWalker();
  void createWalkerInplace(UPtr<MCPWalker>& walker_ptr);
  void allocateWalkerStuffInplace(int walker_index);
  /** }@ */

  void createWalkers();
  void createWalkers(IndexType num_walkers);
  void createWalkers(int num_crowds_,
                     int num_walkers_per_crowd_,
                     IndexType num_walkers,
                     const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& pos);


  /** puts walkers and their "cloned" things into groups in a somewhat general way
   *
   *  Should compile only if ITER is a proper input ITERATOR
   *  Will crash if ITER does point to a std::unique_ptr<WALKER_CONSUMER>
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

  /** The number of cases in which this and get_num_local_walkers is so few that
   *  I strongly suspect it is a design issue.
   *
   *  get_active_walkers is onlyh useful between the setting of num_local_walkers_ and
   *  creation of walkers.  I would reason this is actually a time over which the MCPopulation object
   *  is invalid. Ideally MCPopulation not process any calls in this state, next best would be to only
   *  process calls to become valid.
   */
  //IndexType get_active_walkers() const { return walkers_.size(); }
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
  std::vector<IndexType> get_num_local_walkers_per_node() const { return num_local_walkers_per_node_; }
  void syncWalkersPerNode(Communicate* comm);
  void set_num_global_walkers(IndexType num_global_walkers) { num_global_walkers_ = num_global_walkers; }
  void set_num_local_walkers(IndexType num_local_walkers) { num_local_walkers_ = num_local_walkers; }

  void set_target(IndexType pop) { target_population_ = pop; }
  void set_target_samples(IndexType samples) { target_samples_ = samples; }

  void set_ensemble_property(const MCDataType<QMCTraits::FullPrecRealType>& ensemble_property)
  {
    ensemble_property_ = ensemble_property;
  }

  UPtrVector<MCPWalker>& get_walkers() { return walkers_; }
  UPtrVector<QMCHamiltonian>& get_hamiltonians() { return walker_hamiltonians_; }
  const std::vector<std::pair<int, int>>& get_particle_group_indexes() const { return particle_group_indexes_; }
  const std::vector<RealType>& get_ptclgrp_mass() const { return ptclgrp_mass_; }
  const std::vector<RealType>& get_ptclgrp_inv_mass() const { return ptclgrp_inv_mass_; }
  const std::vector<RealType>& get_ptcl_inv_mass() const { return ptcl_inv_mass_; }

  void set_walker_offsets(std::vector<IndexType> walker_offsets) { walker_offsets_ = walker_offsets; }

  // TODO: the fact this is needed is sad remove need for its existence.
  QMCHamiltonian& get_golden_hamiltonian() { return *hamiltonian_; }
  /** }@ */
};


} // namespace qmcplusplus

#endif
