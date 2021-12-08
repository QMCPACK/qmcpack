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
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "QMCDrivers/WalkerElementsRef.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Utilities/FairDivide.h"

// forward declaration
namespace optimize
{
struct VariableSet;
}

namespace qmcplusplus
{
// forward declaration
class QMCHamiltonian;
class MCPopulation
{
public:
  using MCPWalker          = Walker<QMCTraits, PtclOnLatticeTraits>;
  using WFBuffer           = MCPWalker::WFBuffer_t;
  using RealType           = QMCTraits::RealType;
  using Properties         = MCPWalker::PropertyContainer_t;
  using IndexType          = QMCTraits::IndexType;
  using FullPrecRealType   = QMCTraits::FullPrecRealType;
  using opt_variables_type = optimize::VariableSet;

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

  // By making this a linked list and creating the crowds at the same time we could get first touch.
  UPtrVector<MCPWalker> walkers_;
  UPtrVector<MCPWalker> dead_walkers_;
  std::vector<std::pair<int, int>> particle_group_indexes_;
  std::vector<RealType> ptclgrp_mass_;
  ///1/Mass per species
  std::vector<RealType> ptclgrp_inv_mass_;
  ///1/Mass per particle
  std::vector<RealType> ptcl_inv_mass_;

  // This is necessary MCPopulation is constructed in a simple call scope in QMCDriverFactory from the legacy MCWalkerConfiguration
  // MCPopulation should have QMCMain scope eventually and the driver will just have a reference to it.
  // Then these too can be references.
  TrialWaveFunction* trial_wf_;
  ParticleSet* elec_particle_set_;
  QMCHamiltonian* hamiltonian_;
  WaveFunctionFactory* wf_factory_;
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

  // reference to the captured WalkerConfigurations
  WalkerConfigurations& walker_configs_ref_;

public:
  /** Temporary constructor to deal with MCWalkerConfiguration be the only source of some information
   *  in QMCDriverFactory.
   */
  MCPopulation(int num_ranks,
               int this_rank,
               WalkerConfigurations& mcwc,
               ParticleSet* elecs,
               TrialWaveFunction* trial_wf,
               WaveFunctionFactory* wf_factory,
               QMCHamiltonian* hamiltonian_);

  ~MCPopulation();
  MCPopulation(MCPopulation&) = delete;
  MCPopulation& operator=(MCPopulation&) = delete;
  MCPopulation(MCPopulation&&)           = default;

  /** @ingroup PopulationControl
   *
   *  State Requirement:
   *   * createWalkers must have been called
   *  @{
   */
  WalkerElementsRef spawnWalker();
  void killWalker(MCPWalker&);
  void killLastWalker();
  /** }@ */

  /** Creates walkers with a clone of the golden electron particle set and golden trial wavefunction
   *
   *  \param[in] num_walkers number of living walkers in initial population
   *  \param[in] reserve multiple above that to reserve >=1.0
   */
  void createWalkers(IndexType num_walkers, RealType reserve = 1.0);

  /** distributes walkers and their "cloned" elements to the elements of a vector
   *  of unique_ptr to "walker_consumers". 
   *
   *  a valid "walker_consumer" has a member function of
   *  void addWalker(MCPWalker& walker, ParticleSet& elecs, TrialWaveFunction& twf, QMCHamiltonian& hamiltonian);
   */
  template<typename WTTV>
  void redistributeWalkers(WTTV& walker_consumers)
  {
    // The type returned here is dependent on the integral type that the walker_consumers
    // use to return there size.
    auto walkers_per_crowd = fairDivide(walkers_.size(), walker_consumers.size());

    auto walker_index = 0;
    for (int i = 0; i < walker_consumers.size(); ++i)
    {
      walker_consumers[i]->clearWalkers();
      for (int j = 0; j < walkers_per_crowd[i]; ++j)
      {
        walker_consumers[i]->addWalker(*walkers_[walker_index], *walker_elec_particle_sets_[walker_index],
                                       *walker_trial_wavefunctions_[walker_index], *walker_hamiltonians_[walker_index]);
        ++walker_index;
      }
    }
  }

  void syncWalkersPerRank(Communicate* comm);
  void measureGlobalEnergyVariance(Communicate& comm, FullPrecRealType& ener, FullPrecRealType& variance) const;

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

  // accessor to the gold copy
  const ParticleSet* get_golden_electrons() const { return elec_particle_set_; }
  ParticleSet* get_golden_electrons() { return elec_particle_set_; }
  const TrialWaveFunction& get_golden_twf() const { return *trial_wf_; }
  TrialWaveFunction& get_golden_twf() { return *trial_wf_; }
  // TODO: the fact this is needed is sad remove need for its existence.
  QMCHamiltonian& get_golden_hamiltonian() { return *hamiltonian_; }
  WaveFunctionFactory& get_wf_factory() { return *wf_factory_; }
  
  void set_num_global_walkers(IndexType num_global_walkers) { num_global_walkers_ = num_global_walkers; }
  void set_num_local_walkers(IndexType num_local_walkers) { num_local_walkers_ = num_local_walkers; }

  void set_target(IndexType pop) { target_population_ = pop; }
  void set_target_samples(IndexType samples) { target_samples_ = samples; }

  void set_ensemble_property(const MCDataType<QMCTraits::FullPrecRealType>& ensemble_property)
  {
    ensemble_property_ = ensemble_property;
  }

  UPtrVector<MCPWalker>& get_walkers() { return walkers_; }
  const UPtrVector<MCPWalker>& get_walkers() const { return walkers_; }
  const UPtrVector<MCPWalker>& get_dead_walkers() const { return dead_walkers_; }

  UPtrVector<ParticleSet>& get_elec_particle_sets() { return walker_elec_particle_sets_; }

  UPtrVector<TrialWaveFunction>& get_twfs() { return walker_trial_wavefunctions_; }
  UPtrVector<TrialWaveFunction>& get_dead_twfs() { return dead_walker_trial_wavefunctions_; }

  UPtrVector<QMCHamiltonian>& get_hamiltonians() { return walker_hamiltonians_; }
  UPtrVector<QMCHamiltonian>& get_dead_hamiltonians() { return dead_walker_hamiltonians_; }

  /** Non threadsafe access to walkers and their elements
   *  
   *  Prefer to distribute the walker elements and access
   *  through a crowd to support the concurrency design.
   *
   *  You should not use this unless absolutely necessary.
   *  That doesn't include that you would rather just use
   *  omp parallel and ignore concurrency.
   */
  WalkerElementsRef getWalkerElementsRef(const size_t walker_index);

  /** As long as walker WalkerElements is used we need this for unit tests
   *
   *  As operator[] don't use it to ignore the concurrency design.
   */
  std::vector<WalkerElementsRef> get_walker_elements();

  const std::vector<std::pair<int, int>>& get_particle_group_indexes() const { return particle_group_indexes_; }
  const std::vector<RealType>& get_ptclgrp_mass() const { return ptclgrp_mass_; }
  const std::vector<RealType>& get_ptclgrp_inv_mass() const { return ptclgrp_inv_mass_; }
  const std::vector<RealType>& get_ptcl_inv_mass() const { return ptcl_inv_mass_; }

  /** }@ */


  /// Set variational parameters for the per-walker copies of the wavefunction.
  void set_variational_parameters(const opt_variables_type& active);

  /// check if all the internal vector contain consistent sizes;
  void checkIntegrity() const;

  WalkerConfigurations& getWalkerConfigsRef() { return walker_configs_ref_; }

  // save walker configurations to walker_configs_ref_
  void saveWalkerConfigurations();
};

} // namespace qmcplusplus

#endif
