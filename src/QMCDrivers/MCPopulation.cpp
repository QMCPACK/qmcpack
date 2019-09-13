//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: MCWalkerConfiguration.cpp, QMCUpdate.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/MCPopulation.h"

#include "Configuration.h"
#include "Concurrency/TasksOneToOne.hpp"

namespace qmcplusplus
{
MCPopulation::MCPopulation(MCWalkerConfiguration& mcwc,
                           ParticleSet* elecs,
                           TrialWaveFunction* trial_wf,
                           QMCHamiltonian* hamiltonian)
    : trial_wf_(trial_wf), elec_particle_set_(elecs), hamiltonian_(hamiltonian)
{
  walker_offsets_     = mcwc.WalkerOffsets;
  num_global_walkers_ = mcwc.GlobalNumWalkers;
  num_local_walkers_  = mcwc.LocalNumWalkers;
  num_particles_      = mcwc.getParticleNum();
  size_dataset        = mcwc.DataSet.size();
  
  // MCWalkerConfiguration doesn't give actual number of groups
  num_groups_ = mcwc.groups();
  particle_group_indexes_.resize(num_groups_);
  for (int i = 0; i < num_groups_; ++i)
  {
    particle_group_indexes_[i].first  = mcwc.first(i);
    particle_group_indexes_[i].second = mcwc.last(i);
  }
  ptclgrp_mass_.resize(num_groups_);
  for (int ig = 0; ig < num_groups_; ++ig)
    ptclgrp_mass_[ig] = mcwc.Mass[ig];
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

/** Default creates walkers equal to num_local_walkers_ and zeroed positions
 */
void MCPopulation::createWalkers()
{
  createWalkers(num_local_walkers_, ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>(num_particles_));
}

/** Creates walkers with starting positions pos and a clone of the electron particle set and trial wavefunction
 *  
 *  Needed
 *  in: DataSet.size()
 *  in.
 */
void MCPopulation::createWalkers(IndexType num_walkers,
                                 const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& positions)
{
  walkers_.resize(num_walkers);

  std::for_each(walkers_.begin(), walkers_.end(), [this, positions](std::unique_ptr<MCPWalker>& walker_ptr) {
    walker_ptr = std::make_unique<MCPWalker>(num_particles_);
    walker_ptr->R.resize(num_particles_);
    walker_ptr->R = positions;
    walker_ptr->registerData();
  });

  // Sadly the wfc makeClone interface depends on the full particle set as a way to not to keep track
  // of what different wave function components depend on. I'm going to try and create a hollow elec PS
  // with an eye toward removing the ParticleSet dependency of WFC components in the future.
  walker_elec_particle_sets_.resize(num_walkers);
  std::for_each(walker_elec_particle_sets_.begin(), walker_elec_particle_sets_.end(),
                [this, positions](std::unique_ptr<ParticleSet>& elec_ps_ptr) {
                  elec_ps_ptr.reset(new ParticleSet(*elec_particle_set_));
                  
                  elec_ps_ptr->R = positions;
                });

  auto it_weps = walker_elec_particle_sets_.begin();
  walker_trial_wavefunctions_.resize(num_walkers);
  auto it_wtw = walker_trial_wavefunctions_.begin();
  walker_hamiltonians_.resize(num_walkers);
  auto it_ham = walker_hamiltonians_.begin();
  while (it_wtw != walker_trial_wavefunctions_.end())
  {
    it_wtw->reset(trial_wf_->makeClone(**it_weps));
    it_ham->reset(hamiltonian_->makeClone(**it_weps,**it_wtw));
    ++it_weps;
    ++it_wtw;
    ++it_ham;
  }

  RefVector<WFBuffer> mcp_wfbuffers;
  mcp_wfbuffers.reserve(num_walkers);
  std::for_each(walkers_.begin(),
                walkers_.end(),
                [&mcp_wfbuffers](auto& walker) {
                  mcp_wfbuffers.push_back((**walker).DataSet);
                });
  
  TrialWaveFunction::flex_registerData(walker_trial_wavefunctions_,
                                       walker_elec_particle_sets_,
                                       mcp_wfbuffers);

  std::for_each(walkers_.begin(),
                walkers_.end(),
                [](auto& walker) {
                  (**walker).DataSet.allocate();
                });

  

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

//   TasksOneToOne<> do_per_crowd(num_crowds);

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


} // namespace qmcplusplus
