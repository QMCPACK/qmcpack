#include "QMCHamiltonians/ChiesaCorrection.h"

namespace qmcplusplus {
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
}
