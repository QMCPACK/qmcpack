//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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
#include "Utilities/FairDivide.h"
namespace qmcplusplus
{
class MCPopulation
{
public:
  using MCPWalker        = Walker<QMCTraits, PtclOnLatticeTraits>;
  using WFBuffer         = MCPWalker::WFBuffer_t;
  using RealType         = QMCTraits::RealType;
  using Properties       = MCPWalker::PropertyContainer_t;
  using IndexType        = QMCTraits::IndexType;
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  
private:
  // Potential thread safety issue
  MCDataType<QMCTraits::FullPrecRealType> ensemble_property_;

  IndexType num_global_walkers_ = 0;
  IndexType num_local_walkers_  = 0;
  IndexType num_particles_      = 0;
  IndexType num_groups_         = 0;
  IndexType max_samples_        = 0;
  IndexType target_population_  = 0;
  IndexType target_samples_     = 0;
  //Properties properties_;
  ParticleSet ions_;

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
  // MCPopulation should have QMCMain scope eventually and the driver will just have a reference to it.
  TrialWaveFunction* trial_wf_;
  ParticleSet* elec_particle_set_;
  QMCHamiltonian* hamiltonian_;
  // At the moment these are "clones" but I think this design pattern smells.
  UPtrVector<ParticleSet> walker_elec_particle_sets_;
  UPtrVector<TrialWaveFunction> walker_trial_wavefunctions_;
  UPtrVector<QMCHamiltonian> walker_hamiltonians_;

  // We still haven't cleaned up the dependence between different walker elements so they all need to be tracked
  // as in the legacy implementation.
  UPtrVector<ParticleSet> dead_walker_elec_particle_sets_;
  UPtrVector<TrialWaveFunction> dead_walker_trial_wavefunctions_;
  UPtrVector<QMCHamiltonian> dead_walker_hamiltonians_;

  // MCPopulation immutables
  // would be nice if they were const but we'd lose the default move assignment
  int num_ranks_;
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
               int this_rank);

  MCPopulation(int num_ranks,
               ParticleSet* elecs,
               TrialWaveFunction* trial_wf,
               QMCHamiltonian* hamiltonian,
               int this_rank);

  MCPopulation(MCPopulation&) = delete;
  MCPopulation& operator=(MCPopulation&) = delete;
  MCPopulation(MCPopulation&&)           = default;

  /** @ingroup PopulationControl
   *
   *  State Requirement:
   *   * createWalkers must have been called
   *  @{
   */
  MCPWalker* spawnWalker();
  void killWalker(MCPWalker&);
  void killLastWalker();
  void createWalkerInplace(UPtr<MCPWalker>& walker_ptr);
  void allocateWalkerStuffInplace(int walker_index);
  /** }@ */

  /** Creates walkers with a clone of the golden electron particle set and golden trial wavefunction
   *
   *  \param[in] num_walkers number of living walkers in initial population
   *  \param[in] reserve multiple above that to reserve >=1.0
   */
  void createWalkers(IndexType num_walkers,RealType reserve = 1.0);
  void createWalkers(int num_crowds_,
                     int num_walkers_per_crowd_,
                     IndexType num_walkers,
                     const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& pos);


  /** puts walkers and their "cloned" things into groups in a somewhat general way
   *
   *  Should compile only if ITER is a proper input ITERATOR
   *  Will crash if ITER does point to a std::unqiue_ptr<WALKER_CONSUMER>
   *
   *  The logic here to minimize moves of walkers from one crowd to another
   *  When num_walkers % walkers_per_crowd is true then at the end the extra
   *  walkers are distributed one by one to crowds.
   *
   */
  template<typename WTTV>
  void distributeWalkers(WTTV& walker_consumer)
  {
    auto walkers_per_crowd = fairDivide(walkers_.size(), walker_consumer.size());

    auto walker_index = 0;
    for (int i = 0; i < walker_consumer.size(); ++i)
    {
      for(int j = 0; j < walkers_per_crowd[i]; ++j)
      {
        walker_consumer[i]->addWalker(*walkers_[walker_index], *walker_elec_particle_sets_[walker_index], *walker_trial_wavefunctions_[walker_index], *walker_hamiltonians_[walker_index]);
        ++walker_index;
      }
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
  int get_rank() const { return rank_; }
  IndexType get_num_global_walkers() const { return num_global_walkers_; }
  IndexType get_num_local_walkers() const { return num_local_walkers_; }
  IndexType get_num_particles() const { return num_particles_; }
  IndexType get_max_samples() const { return max_samples_; }
  IndexType get_target_population() const { return target_population_; }
  IndexType get_target_samples() const { return target_samples_; }
  //const Properties& get_properties() const { return properties_; }
  const SpeciesSet& get_species_set() const { return species_set_; }
  const ParticleSet& get_ions() const { return ions_; }
  const ParticleSet* get_golden_electrons() const { return elec_particle_set_; }
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
  UPtrVector<MCPWalker>& get_dead_walkers() { return dead_walkers_; }

  UPtrVector<QMCHamiltonian>& get_hamiltonians() { return walker_hamiltonians_; }
  UPtrVector<QMCHamiltonian>& get_dead_hamiltonians() { return dead_walker_hamiltonians_; }

  const std::vector<std::pair<int, int>>& get_particle_group_indexes() const { return particle_group_indexes_; }
  const std::vector<RealType>& get_ptclgrp_mass() const { return ptclgrp_mass_; }
  const std::vector<RealType>& get_ptclgrp_inv_mass() const { return ptclgrp_inv_mass_; }
  const std::vector<RealType>& get_ptcl_inv_mass() const { return ptcl_inv_mass_; }

  // TODO: the fact this is needed is sad remove need for its existence.
  QMCHamiltonian& get_golden_hamiltonian() { return *hamiltonian_; }
  /** }@ */


  /// Set variational parameters for the per-walker copies of the wavefunction.
  void set_variational_parameters(const opt_variables_type& active);
};

} // namespace qmcplusplus

#endif
