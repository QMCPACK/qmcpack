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
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
MCPopulation::MCPopulation(int num_ranks,
                           int this_rank,
                           WalkerConfigurations& mcwc,
                           ParticleSet* elecs,
                           TrialWaveFunction* trial_wf,
                           QMCHamiltonian* hamiltonian)
    : trial_wf_(trial_wf),
      elec_particle_set_(elecs),
      hamiltonian_(hamiltonian),
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

/** Default creates walkers equal to num_local_walkers_ and zeroed positions
 */
//void MCPopulation::createWalkers() { createWalkers(num_local_walkers_); }

/** we could also search for walker_ptr
 */
void MCPopulation::allocateWalkerStuffInplace(int walker_index)
{
  walker_trial_wavefunctions_[walker_index]->registerData(*(walker_elec_particle_sets_[walker_index]),
                                                          walkers_[walker_index]->DataSet);
  walkers_[walker_index]->DataSet.allocate();
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
  auto createWalker = [this](UPtr<MCPWalker>& walker_ptr) {
    walker_ptr        = std::make_unique<MCPWalker>(num_particles_);
    walker_ptr->R     = elec_particle_set_->R;
    walker_ptr->spins = elec_particle_set_->spins;
    // Side effect of this changes size of walker_ptr->Properties if done after registerData() you end up with
    // a bad buffer.
    walker_ptr->Properties = elec_particle_set_->Properties;
    walker_ptr->registerData();
  };

  for (auto& walker_ptr : walkers_)
    createWalker(walker_ptr);

  for (int iw = 0; iw < std::min(walkers_.size(), walker_configs_ref_.WalkerList.size()); iw++)
    *walkers_[iw] = *walker_configs_ref_[iw];

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

  outputManager.pause();

  // Sadly the wfc makeClone interface depends on the full particle set as a way to not to keep track
  // of what different wave function components depend on. I'm going to try and create a hollow elec PS
  // with an eye toward removing the ParticleSet dependency of WFC components in the future.
  walker_elec_particle_sets_.resize(num_walkers_plus_reserve);
  std::for_each(walker_elec_particle_sets_.begin(), walker_elec_particle_sets_.end(),
                [this](std::unique_ptr<ParticleSet>& elec_ps_ptr) {
                  elec_ps_ptr = std::make_unique<ParticleSet>(*elec_particle_set_);
                });

  auto it_weps = walker_elec_particle_sets_.begin();
  walker_trial_wavefunctions_.resize(num_walkers_plus_reserve);
  auto it_wtw = walker_trial_wavefunctions_.begin();
  walker_hamiltonians_.resize(num_walkers_plus_reserve);
  auto it_ham = walker_hamiltonians_.begin();
  while (it_wtw != walker_trial_wavefunctions_.end())
  {
    it_wtw->reset(trial_wf_->makeClone(**it_weps));
    it_ham->reset(hamiltonian_->makeClone(**it_weps, **it_wtw));
    ++it_weps;
    ++it_wtw;
    ++it_ham;
  }

  outputManager.resume();

  RefVector<WFBuffer> mcp_wfbuffers;
  mcp_wfbuffers.reserve(num_walkers_plus_reserve);
  std::for_each(walkers_.begin(), walkers_.end(),
                [&mcp_wfbuffers](auto& walker) { mcp_wfbuffers.push_back((*walker).DataSet); });

  TrialWaveFunction::flex_registerData(walker_trial_wavefunctions_, walker_elec_particle_sets_, mcp_wfbuffers);

  std::for_each(walkers_.begin(), walkers_.end(), [](auto& walker) {
    MCPWalker& this_walker = *walker;
    this_walker.DataSet.allocate();
  });

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
    MCPWalker last_walker = *(walkers_.back());
    walkers_.push_back(std::make_unique<MCPWalker>(last_walker));

    // There is no value in doing this here because its going to be wiped out
    // When we load from the receive buffer. It also won't necessarily be correct
    // Because the buffer is changed by Hamiltonians and wavefunctions that
    // Add to the dataSet.

    //walkers_.back()->R          = elec_particle_set_->R;
    //walkers_.back()->Properties = elec_particle_set_->Properties;
    //walkers_.back()->registerData();

    walker_elec_particle_sets_.emplace_back(std::make_unique<ParticleSet>(*elec_particle_set_));
    walker_trial_wavefunctions_.push_back(UPtr<TrialWaveFunction>{});
    walker_trial_wavefunctions_.back().reset(trial_wf_->makeClone(*(walker_elec_particle_sets_.back())));
    walker_hamiltonians_.push_back(UPtr<QMCHamiltonian>{});
    walker_hamiltonians_.back().reset(
        hamiltonian_->makeClone(*(walker_elec_particle_sets_.back()), *(walker_trial_wavefunctions_.back())));
    // Dito
    //walker_trial_wavefunctions_.back()->registerData(*(walker_elec_particle_sets_.back()), walkers_.back()->DataSet);
    //walkers_.back()->DataSet.allocate();
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
      (*it_walkers)->DataSet.zero();
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

void MCPopulation::syncWalkersPerNode(Communicate* comm)
{
  std::vector<IndexType> num_local_walkers_per_node(comm->size(), 0);
  ;

  num_local_walkers_per_node[comm->rank()] = num_local_walkers_;
  comm->allreduce(num_local_walkers_per_node);

  num_global_walkers_ = std::accumulate(num_local_walkers_per_node.begin(), num_local_walkers_per_node.end(), 0);
}


void MCPopulation::set_variational_parameters(const opt_variables_type& active)
{
  for (auto it_twfs = walker_trial_wavefunctions_.begin(); it_twfs != walker_trial_wavefunctions_.end(); ++it_twfs)
  {
    (*it_twfs).get()->resetParameters(active);
  }
}

/** Creates walkers doing their first touch in their crowd (thread) context
 *
 *  This is basically premature optimization but I wanted to check if this sort of thing
 *  would work.  It seems to. This sort of structure not an #omp parallel section must be used in the driver.
 *  No new bare openmp directives should be added to code with updated design.
 */
// void MCPopulation::createWalkers(int num_crowds,
//                                  int walkers_per_crowd,
//                                  IndexType num_walkers,
//                                  const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& positions)
// {
//   walkers_.resize(num_walkers);

//   ParallelExecutor<> do_per_crowd(num_crowds);

//   std::vector<std::unique_ptr<std::vector<std::unique_ptr<MCPWalker>>>> walkers_per_crowd_per_slot;
//   walkers_per_crowd_per_slot.resize(num_crowds);
//   auto first_touch_create_walkers =
//       [this, &walkers_per_crowd,
//        &positions](int crowd_index, std::vector<std::unique_ptr<std::vector<std::unique_ptr<MCPWalker>>>>& wpcps) {
//         wpcps[crowd_index] = std::make_unique<std::vector<std::unique_ptr<MCPWalker>>>(walkers_per_crowd);
//         std::vector<std::unique_ptr<MCPWalker>>& this_crowds_walkers = *(wpcps[crowd_index]);
//         this_crowds_walkers.resize(walkers_per_crowd);
//         for (int i = 0; i < walkers_per_crowd; ++i)
//         {
//           std::unique_ptr<MCPWalker>& walker_uptr = this_crowds_walkers[i];
//           walker_uptr.reset(new MCPWalker(num_particles_));
//           walker_uptr->R.resize(num_particles_);
//           walker_uptr->R = positions;
//         }
//       };
//   do_per_crowd(first_touch_create_walkers, walkers_per_crowd_per_slot);

//   auto walkers_it = walkers_.begin();
//   std::for_each(walkers_per_crowd_per_slot.begin(), walkers_per_crowd_per_slot.end(),
//                 [&walkers_it](std::unique_ptr<std::vector<std::unique_ptr<MCPWalker>>>& per_crowd_ptr) {
//                   std::vector<std::unique_ptr<MCPWalker>>& walkers_per_crowd = *per_crowd_ptr;
//                   for (int i = 0; i < walkers_per_crowd.size(); ++i)
//                   {
//                     *walkers_it = std::move(walkers_per_crowd[i]);
//                     ++walkers_it;
//                   }
//                 });
// }

void MCPopulation::saveWalkerConfigurations()
{
  walker_configs_ref_.resize(walker_elec_particle_sets_.size(), elec_particle_set_->getTotalNum());
  for (int iw = 0; iw < walker_elec_particle_sets_.size(); iw++)
    walker_elec_particle_sets_[iw]->saveWalker(*walker_configs_ref_[iw]);
}


} // namespace qmcplusplus
