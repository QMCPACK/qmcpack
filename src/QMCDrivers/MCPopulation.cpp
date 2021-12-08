//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: MCWalkerConfiguration.cpp, QMCUpdate.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <numeric>

#include "MCPopulation.h"
#include "Configuration.h"
#include "Concurrency/ParallelExecutor.hpp"
#include "Message/CommOperators.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
MCPopulation::MCPopulation(int num_ranks,
                           int this_rank,
                           WalkerConfigurations& mcwc,
                           ParticleSet* elecs,
                           TrialWaveFunction* trial_wf,
                           WaveFunctionFactory* wf_factory,
                           QMCHamiltonian* hamiltonian)
    : trial_wf_(trial_wf),
      elec_particle_set_(elecs),
      hamiltonian_(hamiltonian),
      wf_factory_(wf_factory),
      num_ranks_(num_ranks),
      rank_(this_rank),
      walker_configs_ref_(mcwc)
{
  num_global_walkers_ = mcwc.getGlobalNumWalkers();
  num_local_walkers_  = mcwc.getActiveWalkers();
  num_particles_      = elecs->getTotalNum();

  // MCWalkerConfiguration doesn't give actual number of groups
  num_groups_ = elecs->groups();
  particle_group_indexes_.resize(num_groups_);
  for (int i = 0; i < num_groups_; ++i)
  {
    particle_group_indexes_[i].first  = elecs->first(i);
    particle_group_indexes_[i].second = elecs->last(i);
  }
  ptclgrp_mass_.resize(num_groups_);
  for (int ig = 0; ig < num_groups_; ++ig)
    ptclgrp_mass_[ig] = elecs->Mass[ig];
  ptclgrp_inv_mass_.resize(num_groups_);
  for (int ig = 0; ig < num_groups_; ++ig)
    ptclgrp_inv_mass_[ig] = 1.0 / ptclgrp_mass_[ig];
  ptcl_inv_mass_.resize(num_particles_);
  for (int ig = 0; ig < num_groups_; ++ig)
  {
    for (int iat = particle_group_indexes_[ig].first; iat < particle_group_indexes_[ig].second; ++iat)
      ptcl_inv_mass_[iat] = ptclgrp_inv_mass_[ig];
  }
}

MCPopulation::~MCPopulation()
{
  // if there are active walkers, save them to lightweight walker configuration list.
  if (walkers_.size())
    saveWalkerConfigurations();
}

void MCPopulation::createWalkers(IndexType num_walkers, RealType reserve)
{
  IndexType num_walkers_plus_reserve = static_cast<IndexType>(num_walkers * reserve);

  // Hack to hopefully insure no truly new walkers will be made by spawn, since I suspect that
  // doesn't capture everything that needs to make a walker + elements valid to load from a transferred
  // buffer;
  // Ye: need to resize walker_t and ParticleSet Properties
  // Really MCPopulation does not own this elec_particle_set_  seems like it should be immutable
  elec_particle_set_->Properties.resize(1, elec_particle_set_->PropertyList.size());

  // This pattern is begging for a micro benchmark, is this really better
  // than the simpler walkers_.pushback;
  walkers_.resize(num_walkers_plus_reserve);
  walker_elec_particle_sets_.resize(num_walkers_plus_reserve);
  walker_trial_wavefunctions_.resize(num_walkers_plus_reserve);
  walker_hamiltonians_.resize(num_walkers_plus_reserve);

  outputManager.pause();

  //this part is time consuming, it must be threaded and calls should be thread-safe.
#pragma omp parallel for
  for (size_t iw = 0; iw < num_walkers_plus_reserve; iw++)
  {
    walkers_[iw]             = std::make_unique<MCPWalker>(num_particles_);
    walkers_[iw]->R          = elec_particle_set_->R;
    walkers_[iw]->spins      = elec_particle_set_->spins;
    walkers_[iw]->Properties = elec_particle_set_->Properties;
    walkers_[iw]->registerData();
    walkers_[iw]->DataSet.allocate();

    if (iw < walker_configs_ref_.WalkerList.size())
      *walkers_[iw] = *walker_configs_ref_[iw];

    walker_elec_particle_sets_[iw]  = std::make_unique<ParticleSet>(*elec_particle_set_);
    walker_trial_wavefunctions_[iw] = trial_wf_->makeClone(*walker_elec_particle_sets_[iw]);
    walker_hamiltonians_[iw] =
        hamiltonian_->makeClone(*walker_elec_particle_sets_[iw], *walker_trial_wavefunctions_[iw]);
  };

  outputManager.resume();

  int num_walkers_created = 0;
  for (auto& walker_ptr : walkers_)
  {
    if (walker_ptr->ID == 0)
    {
      // And so walker ID's start at one because 0 is magic.
      // \todo This is C++ all indexes start at 0, make uninitialized ID = -1
      walker_ptr->ID       = (++num_walkers_created) * num_ranks_ + rank_;
      walker_ptr->ParentID = walker_ptr->ID;
    }
  }

  // kill and spawn walkers update the state variable num_local_walkers_
  // so it must start at the number of reserved walkers
  num_local_walkers_ = num_walkers_plus_reserve;

  IndexType extra_walkers = num_walkers_plus_reserve - num_walkers;
  // Now we kill the extra reserve walkers and elements that we made.
  for (int i = 0; i < extra_walkers; ++i)
    killLastWalker();
}

WalkerElementsRef MCPopulation::getWalkerElementsRef(const size_t index)
{
  return {*walkers_[index], *walker_elec_particle_sets_[index], *walker_trial_wavefunctions_[index]};
}

std::vector<WalkerElementsRef> MCPopulation::get_walker_elements()
{
  std::vector<WalkerElementsRef> walker_elements;
  for (int iw = 0; iw < walkers_.size(); ++iw)
  {
    walker_elements.emplace_back(*walkers_[iw], *walker_elec_particle_sets_[iw], *walker_trial_wavefunctions_[iw]);
  }
  return walker_elements;
}

/** creates a walker and returns a reference
 *
 *  Walkers are reused
 *  It would be better if this could be done just by
 *  reusing memory.
 *  Not thread safe.
 */
WalkerElementsRef MCPopulation::spawnWalker()
{
  ++num_local_walkers_;
  outputManager.pause();

  if (dead_walkers_.size() > 0)
  {
    walkers_.push_back(std::move(dead_walkers_.back()));
    dead_walkers_.pop_back();
    walker_elec_particle_sets_.push_back(std::move(dead_walker_elec_particle_sets_.back()));
    dead_walker_elec_particle_sets_.pop_back();
    walker_trial_wavefunctions_.push_back(std::move(dead_walker_trial_wavefunctions_.back()));
    dead_walker_trial_wavefunctions_.pop_back();
    walker_hamiltonians_.push_back(std::move(dead_walker_hamiltonians_.back()));
    dead_walker_hamiltonians_.pop_back();
    // Emulating the legacy implementation valid walker elements were created with the initial walker and DataSet
    // registration and allocation were done then so are not necessary when resurrecting walkers and elements
    walkers_.back()->Generation         = 0;
    walkers_.back()->Age                = 0;
    walkers_.back()->ReleasedNodeWeight = 1.0;
    walkers_.back()->ReleasedNodeAge    = 0;
    walkers_.back()->Multiplicity       = 1.0;
    walkers_.back()->Weight             = 1.0;
  }
  else
  {
    app_warning() << "Spawning walker number " << walkers_.size() + 1
                  << " outside of reserves, this ideally should never happend." << std::endl;
    walkers_.push_back(std::make_unique<MCPWalker>(*(walkers_.back())));

    // There is no value in doing this here because its going to be wiped out
    // When we load from the receive buffer. It also won't necessarily be correct
    // Because the buffer is changed by Hamiltonians and wavefunctions that
    // Add to the dataSet.

    walker_elec_particle_sets_.emplace_back(std::make_unique<ParticleSet>(*elec_particle_set_));
    walker_trial_wavefunctions_.emplace_back(trial_wf_->makeClone(*walker_elec_particle_sets_.back()));
    walker_hamiltonians_.emplace_back(
        hamiltonian_->makeClone(*walker_elec_particle_sets_.back(), *walker_trial_wavefunctions_.back()));
    walkers_.back()->Multiplicity = 1.0;
    walkers_.back()->Weight       = 1.0;
  }

  outputManager.resume();
  return {*walkers_.back().get(), *walker_elec_particle_sets_.back().get(), *walker_trial_wavefunctions_.back().get()};
}

/** Kill last walker (just barely)
 *
 *  By kill we mean put it and all its elements in a "dead" list.
 */
void MCPopulation::killLastWalker()
{
  --num_local_walkers_;
  // kill the walker but just barely we need all its setup and connections to remain
  dead_walkers_.push_back(std::move(walkers_.back()));
  walkers_.pop_back();
  dead_walker_elec_particle_sets_.push_back(std::move(walker_elec_particle_sets_.back()));
  walker_elec_particle_sets_.pop_back();
  dead_walker_trial_wavefunctions_.push_back(std::move(walker_trial_wavefunctions_.back()));
  walker_trial_wavefunctions_.pop_back();
  dead_walker_hamiltonians_.push_back(std::move(walker_hamiltonians_.back()));
  walker_hamiltonians_.pop_back();
}
/** Kill a walker (just barely)
 *
 *  By kill we mean put it and all its elements in a "dead" list.
 */
void MCPopulation::killWalker(MCPWalker& walker)
{
  // find the walker and move its pointer to the dead walkers vector
  auto it_walkers = walkers_.begin();
  auto it_psets   = walker_elec_particle_sets_.begin();
  auto it_twfs    = walker_trial_wavefunctions_.begin();
  auto it_hams    = walker_hamiltonians_.begin();
  while (it_walkers != walkers_.end())
  {
    if (&walker == (*it_walkers).get())
    {
      dead_walkers_.push_back(std::move(*it_walkers));
      walkers_.erase(it_walkers);
      dead_walker_elec_particle_sets_.push_back(std::move(*it_psets));
      walker_elec_particle_sets_.erase(it_psets);
      dead_walker_trial_wavefunctions_.push_back(std::move(*it_twfs));
      walker_trial_wavefunctions_.erase(it_twfs);
      dead_walker_hamiltonians_.push_back(std::move(*it_hams));
      walker_hamiltonians_.erase(it_hams);
      --num_local_walkers_;
      return;
    }
    ++it_walkers;
    ++it_psets;
    ++it_twfs;
    ++it_hams;
  }
  throw std::runtime_error("Attempt to kill nonexistent walker in MCPopulation!");
}

void MCPopulation::syncWalkersPerRank(Communicate* comm)
{
  std::vector<IndexType> num_local_walkers_per_rank(comm->size(), 0);

  num_local_walkers_per_rank[comm->rank()] = num_local_walkers_;
  comm->allreduce(num_local_walkers_per_rank);

  num_global_walkers_ = std::accumulate(num_local_walkers_per_rank.begin(), num_local_walkers_per_rank.end(), 0);
}

void MCPopulation::measureGlobalEnergyVariance(Communicate& comm,
                                               FullPrecRealType& ener,
                                               FullPrecRealType& variance) const
{
  std::vector<FullPrecRealType> weight_energy_variance(3, 0.0);
  for (int iw = 0; iw < walker_elec_particle_sets_.size(); iw++)
  {
    auto w = walkers_[iw]->Weight;
    auto e = walker_hamiltonians_[iw]->getLocalEnergy();
    weight_energy_variance[0] += w;
    weight_energy_variance[1] += w * e;
    weight_energy_variance[2] += w * e * e;
  }

  comm.allreduce(weight_energy_variance);
  ener     = weight_energy_variance[1] / weight_energy_variance[0];
  variance = weight_energy_variance[2] / weight_energy_variance[0] - ener * ener;
}

void MCPopulation::set_variational_parameters(const opt_variables_type& active)
{
  for (auto it_twfs = walker_trial_wavefunctions_.begin(); it_twfs != walker_trial_wavefunctions_.end(); ++it_twfs)
  {
    (*it_twfs).get()->resetParameters(active);
  }
}

void MCPopulation::checkIntegrity() const
{
  // check active walkers
  const size_t num_local_walkers_active = num_local_walkers_;
  if (walkers_.size() != num_local_walkers_active)
    throw std::runtime_error("walkers_ has inconsistent size");
  if (walker_elec_particle_sets_.size() != num_local_walkers_active)
    throw std::runtime_error("walker_elec_particle_sets_ has inconsistent size");
  if (walker_trial_wavefunctions_.size() != num_local_walkers_active)
    throw std::runtime_error("walker_trial_wavefunctions_ has inconsistent size");
  if (walker_trial_wavefunctions_.size() != num_local_walkers_active)
    throw std::runtime_error("walker_trial_wavefunctions_ has inconsistent size");
  if (walker_hamiltonians_.size() != num_local_walkers_active)
    throw std::runtime_error("walker_hamiltonians_ has inconsistent size");

  // check dead walkers
  const size_t num_local_walkers_dead = dead_walkers_.size();
  if (dead_walker_elec_particle_sets_.size() != num_local_walkers_dead)
    throw std::runtime_error("dead_walker_elec_particle_sets_ has inconsistent size");
  if (dead_walker_trial_wavefunctions_.size() != num_local_walkers_dead)
    throw std::runtime_error("dead_walker_trial_wavefunctions_ has inconsistent size");
  if (dead_walker_trial_wavefunctions_.size() != num_local_walkers_dead)
    throw std::runtime_error("dead_walker_trial_wavefunctions_ has inconsistent size");
  if (dead_walker_hamiltonians_.size() != num_local_walkers_dead)
    throw std::runtime_error("dead_walker_hamiltonians_ has inconsistent size");
}

void MCPopulation::saveWalkerConfigurations()
{
  walker_configs_ref_.resize(walker_elec_particle_sets_.size(), elec_particle_set_->getTotalNum());
  for (int iw = 0; iw < walker_elec_particle_sets_.size(); iw++)
    walker_elec_particle_sets_[iw]->saveWalker(*walker_configs_ref_[iw]);
}


} // namespace qmcplusplus
