#include "Particle/MCPopulation.h"
namespace qmcplusplus
{
MCPopulation::MCPopulation(MCWalkerConfiguration& mcwc)
{
  walker_offsets_     = mcwc.WalkerOffsets;
  num_global_walkers_ = mcwc.GlobalNumWalkers;
  num_local_walkers_  = mcwc.LocalNumWalkers;
}
} // namespace qmcplusplus
