#ifndef QMCPLUSPLUS_RMC_FACTORY_H
#define QMCPLUSPLUS_RMC_FACTORY_H
#include "QMCDrivers/QMCDriver.h"

namespace qmcplusplus
{
  class ParticleSetPool;
  class HamiltonianPool;

  struct RMCFactory
  {
    int RMCMode;
    xmlNodePtr myNode;
      RMCFactory (int vmode, xmlNodePtr cur):RMCMode (vmode), myNode (cur)
    {
    }
    QMCDriver *create (MCWalkerConfiguration & w, TrialWaveFunction & psi,
		       QMCHamiltonian & h, ParticleSetPool & ptclpool,
		       HamiltonianPool & hpool, WaveFunctionPool & ppool);
  };
}

#endif
