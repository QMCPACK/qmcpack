#include "Particle/MCPopulation.h"

#include "Configuration.h"
#include "Concurrency/TasksOneToOne.hpp"

namespace qmcplusplus
{
MCPopulation::MCPopulation(MCWalkerConfiguration& mcwc)
{
  walker_offsets_     = mcwc.WalkerOffsets;
  num_global_walkers_ = mcwc.GlobalNumWalkers;
  num_local_walkers_  = mcwc.LocalNumWalkers;
  num_particles_      = mcwc.getParticleNum();
  // MCWalkerConfiguration doesn't give actual number of groups
  num_groups_         = mcwc.groups() + 1;
  particle_group_indexes_.resize(num_groups_);
  for(int i = 0; i < num_groups_; ++i)
  {
    particle_group_indexes_[i].first = mcwc.first(i);
    particle_group_indexes_[i].second = mcwc.last(i);
  }
}

/** Default creates walkers equal to num_local_walkers_ and zeroed positions
 */
void MCPopulation::createWalkers()
{
  createWalkers(num_local_walkers_, ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>(num_particles_));
}

/** Creates walkers with starting positions pos
 *
 *  Eventually MCPopulation should not depend on ParticleAttrib
 */
void MCPopulation::createWalkers(IndexType num_walkers,
                                 const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& positions)
{
  walkers_.resize(num_walkers);

  std::for_each(walkers_.begin(), walkers_.end(), [this, positions](std::unique_ptr<MCPWalker>& walker_ptr) {
    walker_ptr = std::make_unique<MCPWalker>(num_particles_);
    walker_ptr->R.resize(num_particles_);
    walker_ptr->R = positions;
  });
}

/** Creates walkers with doing there first touch in their crowd (thread) context
 *
 *  This is basically premature optimization but I wanted to check if this sort of thing
 *  would work.  It seems to.
 */
void MCPopulation::createWalkers(int num_crowds,
                                 int walkers_per_crowd,
                                 IndexType num_walkers,
                                 const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& positions)
{
  walkers_.resize(num_walkers);

  TasksOneToOne<> do_per_crowd(num_crowds);

  std::vector<std::unique_ptr<std::vector<std::unique_ptr<MCPWalker>>>> walkers_per_crowd_per_slot;
  walkers_per_crowd_per_slot.resize(num_crowds);
  auto first_touch_create_walkers =
      [this, &walkers_per_crowd,
       &positions](int crowd_index, std::vector<std::unique_ptr<std::vector<std::unique_ptr<MCPWalker>>>>& wpcps) {
        wpcps[crowd_index] = std::make_unique<std::vector<std::unique_ptr<MCPWalker>>>(walkers_per_crowd);
        std::vector<std::unique_ptr<MCPWalker>>& this_crowds_walkers = *(wpcps[crowd_index]);
        this_crowds_walkers.resize(walkers_per_crowd);
        for (int i = 0; i < walkers_per_crowd; ++i)
        {
          std::unique_ptr<MCPWalker>& walker_uptr = this_crowds_walkers[i];
          walker_uptr.reset(new MCPWalker(num_particles_));
          walker_uptr->R.resize(num_particles_);
          walker_uptr->R = positions;
        }
      };
  do_per_crowd(first_touch_create_walkers, walkers_per_crowd_per_slot);

  auto walkers_it = walkers_.begin();
  std::for_each(walkers_per_crowd_per_slot.begin(), walkers_per_crowd_per_slot.end(),
                [&walkers_it](std::unique_ptr<std::vector<std::unique_ptr<MCPWalker>>>& per_crowd_ptr) {
                  std::vector<std::unique_ptr<MCPWalker>>& walkers_per_crowd = *per_crowd_ptr;
                  for (int i = 0; i < walkers_per_crowd.size(); ++i)
                  {
                    *walkers_it = std::move(walkers_per_crowd[i]);
                    ++walkers_it;
                  }
                });
}



} // namespace qmcplusplus
