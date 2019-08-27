#include "Particle/MCPopulation.h"

#include "Configuration.h"
namespace qmcplusplus
{
MCPopulation::MCPopulation(MCWalkerConfiguration& mcwc)
{
  walker_offsets_     = mcwc.WalkerOffsets;
  num_global_walkers_ = mcwc.GlobalNumWalkers;
  num_local_walkers_  = mcwc.LocalNumWalkers;
  num_particles_ = mcwc.getParticleNum();
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
void MCPopulation::createWalkers(IndexType num_walkers, const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& positions)
{
  walkers_.resize(num_walkers);
  std::for_each(walkers_.begin(),
                walkers_.end(),
                [this, positions](std::unique_ptr<MCPWalker>& walker_ptr)
                    {
                      walker_ptr = std::make_unique<MCPWalker>(num_particles_);
                      walker_ptr->R.resize(num_particles_);
                      walker_ptr->R = positions;
                    });
}

} // namespace qmcplusplus
