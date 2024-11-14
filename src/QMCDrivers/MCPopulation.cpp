//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
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
                           ParticleSet* elecs,
                           TrialWaveFunction* trial_wf,
                           QMCHamiltonian* hamiltonian)
    : trial_wf_(trial_wf), elec_particle_set_(elecs), hamiltonian_(hamiltonian), num_ranks_(num_ranks), rank_(this_rank)
{
  const auto num_groups = elecs->groups();
  ptclgrp_mass_.resize(num_groups);
  ptclgrp_inv_mass_.resize(num_groups);
  for (int ig = 0; ig < num_groups; ++ig)
  {
    ptclgrp_mass_[ig]     = elecs->Mass[elecs->first(ig)];
    ptclgrp_inv_mass_[ig] = 1.0 / ptclgrp_mass_[ig];
  }

  ptcl_inv_mass_.resize(elecs->getTotalNum());
  for (int ig = 0; ig < num_groups; ++ig)
    for (int iat = elecs->first(ig); iat < elecs->last(ig); ++iat)
      ptcl_inv_mass_[iat] = ptclgrp_inv_mass_[ig];
}

MCPopulation::~MCPopulation() = default;

void MCPopulation::fissionHighMultiplicityWalkers()
{
  // we need to do this because spawnWalker changes walkers_ so we
  // can't just iterate on that collection.
  auto good_walkers = convertUPtrToRefVector(walkers_);
  for (MCPWalker& good_walker : good_walkers)
  {
    int num_copies = static_cast<int>(good_walker.Multiplicity);
    while (num_copies > 1)
    {
      auto walker_elements = spawnWalker();
      // In the batched version walker ids are set when walkers are born and are unique,
      // parent ids are set to the walker that provided the initial configuration.
      // If the amplified walker was a transfer from another rank its copys get its ID
      // as their parent, Not the walker id the received walker had on its original rank.
      // So walkers don't have to maintain state that they were transfers.
      // We don't need branching here for transfered and nontransfered high multiplicity
      // walkers. The walker assignment operator could avoid writing to the walker_id of
      // left side walker but perhaps the assignment operator is surprising enough as is.
      auto walker_id         = walker_elements.walker.getWalkerID();
      walker_elements.walker = good_walker;
      // copy the copied from walkers id to parent id.
      walker_elements.walker.setParentID(walker_elements.walker.getWalkerID());
      // put the walkers actual id back.
      walker_elements.walker.setWalkerID(walker_id);
      // fix the multiplicity of the new walker
      walker_elements.walker.Multiplicity = 1.0;
      // keep good walker valid.
      good_walker.Multiplicity -= 1.0;
      num_copies--;
    }
  }
}

void MCPopulation::createWalkers(IndexType num_walkers, const WalkerConfigurations& walker_configs, RealType reserve)
{
  assert(reserve >= 1.0);
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

  // nextWalkerID is not thread safe so we need to get our walker id's here
  // before we enter a parallel section.
  std::vector<long> walker_ids(num_walkers_plus_reserve, 0);
  // notice we only get walker_ids for the number of walkers in the walker_configs
  for (size_t iw = 0; iw < num_walkers; iw++)
    walker_ids[iw] = nextWalkerID();

  // this part is time consuming, it must be threaded and calls should be thread-safe.
  // It would make more sense if it was over crowd threads as the thread locality of the walkers
  // would at least initially be "optimal" Depending on the number of OMP threads and implementation
  // this may be equivalent.
#pragma omp parallel for shared(walker_ids)
  for (size_t iw = 0; iw < num_walkers_plus_reserve; iw++)
  {
    // initialize walkers from existing walker_configs
    if (const auto num_existing_walkers = walker_configs.getActiveWalkers())
    {
      walkers_[iw] = std::make_unique<MCPWalker>(*walker_configs[iw % num_existing_walkers]);
      // An outside section context parent ID is multiplied by -1.
      walkers_[iw]->setParentID(-1 * walker_configs[iw % num_existing_walkers]->getWalkerID());
      walkers_[iw]->setWalkerID(walker_ids[iw]);
    }
    else // these are fresh walkers no incoming walkers
    {
      // These walkers are orphans they don't get their intial configuration from a walkerconfig
      // but from the golden particle set.  They get an walker ID of 0;
      walkers_[iw] = std::make_unique<MCPWalker>(walker_ids[iw], 0 /* parent_id */, elec_particle_set_->getTotalNum());
      // Should these get a randomize from source?
      // This seems to be what happens in legacy but its surprisingly opaque there
      // How is it not undesirable to have all these walkers start from the same positions
      walkers_[iw]->R     = elec_particle_set_->R;
      walkers_[iw]->spins = elec_particle_set_->spins;
    }

    walkers_[iw]->Properties = elec_particle_set_->Properties;
    walkers_[iw]->registerData();
    walkers_[iw]->DataSet.allocate();

    walker_elec_particle_sets_[iw]  = std::make_unique<ParticleSet>(*elec_particle_set_);
    walker_trial_wavefunctions_[iw] = trial_wf_->makeClone(*walker_elec_particle_sets_[iw]);
    walker_hamiltonians_[iw] =
        hamiltonian_->makeClone(*walker_elec_particle_sets_[iw], *walker_trial_wavefunctions_[iw]);
  };

  outputManager.resume();

  // kill and spawn walkers update the state variable num_local_walkers_
  // so it must start at the number of reserved walkers
  num_local_walkers_ = num_walkers_plus_reserve;

  IndexType extra_walkers = num_walkers_plus_reserve - num_walkers;
  // Now we kill the extra reserve walkers and elements that we made.
  for (int i = 0; i < extra_walkers; ++i)
    killLastWalker();
  // And now num_local_walkers_ will be correct.
}

long MCPopulation::nextWalkerID() { return num_walkers_created_++ * num_ranks_ + rank_ + 1; }

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
 *  Not thread safe.
 *
 *  In most cases Walkers that are reused have either died during DMC
 *  or were created dead as a reserve. This is due to the dynamic allocation of the Walker
 *  and walker elements being expensive in time.  If still true this is an important optimization.
 *  I think its entirely possible that the walker elements are bloated and if they only included
 *  necessary per walker mutable elements that this entire complication could be removed.
 *
 *  Walker ID's are handed out per independent trajectory (life) not per allocation.
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
    walkers_.back()->Generation   = 0;
    walkers_.back()->Age          = 0;
    walkers_.back()->Multiplicity = 1.0;
    walkers_.back()->Weight       = 1.0;
    // this does count as a walker creation so it gets a new walker id
    walkers_.back()->setWalkerID(nextWalkerID());
  }
  else
  {
    outputManager.resume();
    auto walker_id = nextWalkerID();
    app_debug() << "Spawning a walker (ID " << walker_id << ") triggers living walker number " << walkers_.size()
                << " allocation. This happens when population starts to fluctuate at the begining of a simulation "
                << "but infrequently when the fluctuation stablizes." << std::endl;
    walkers_.push_back(std::make_unique<MCPWalker>(*(walkers_.back()), walker_id, 0));

    outputManager.pause();

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

void MCPopulation::saveWalkerConfigurations(WalkerConfigurations& walker_configs)
{
  walker_configs.resize(walker_elec_particle_sets_.size(), elec_particle_set_->getTotalNum());
  for (int iw = 0; iw < walker_elec_particle_sets_.size(); iw++)
  {
    walker_configs[iw]->R      = walkers_[iw]->R;
    walker_configs[iw]->spins  = walkers_[iw]->spins;
    walker_configs[iw]->G      = walkers_[iw]->G;
    walker_configs[iw]->L      = walkers_[iw]->L;
    walker_configs[iw]->Weight = walkers_[iw]->Weight;
    walker_configs[iw]->setWalkerID(walkers_[iw]->getWalkerID());
    walker_configs[iw]->setParentID(walkers_[iw]->getParentID());
  }
}
} // namespace qmcplusplus
