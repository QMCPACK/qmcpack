#include "Particle/MCPopulation.h"
namespace qmcplusplus
{
  MCPopulation::MCPopulation(MCWalkerConfiguration& mcwc)
  {
    walker_offsets_ = mcwc.WalkerOffsets;
    global_num_walkers_ = mcwc.GlobalNumWalkers;
    local_num_walkers_ = mcwc.LocalNumWalkers;
  }
}
