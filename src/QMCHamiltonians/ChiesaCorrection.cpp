#include "QMCHamiltonians/ChiesaCorrection.h"

namespace qmcplusplus
{
void
ChiesaCorrection::resetTargetParticleSet(ParticleSet& P)
{
}

bool
ChiesaCorrection::put(xmlNodePtr cur)
{
  return true;
}


QMCHamiltonianBase*
ChiesaCorrection::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return new ChiesaCorrection (qp, psi);
}

ChiesaCorrection::Return_t
ChiesaCorrection::evaluate(ParticleSet& P)
{
  return Value = psi_ref.KECorrection();
}
#ifdef QMC_CUDA
void
ChiesaCorrection::addEnergy(MCWalkerConfiguration &W,
                            vector<RealType> &LocalEnergy)
{
  RealType corr = psi_ref.KECorrection();
  vector<Walker_t*> &walkers = W.WalkerList;
  for (int iw=0; iw<LocalEnergy.size(); iw++)
  {
    walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] = corr;
    LocalEnergy[iw] += corr;
  }
}
#endif
}
