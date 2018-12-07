
#ifndef QMCPLUSPLUS_DIFF_COUNTING_JASTROW_ORBITAL_H
#define QMCPLUSPLUS_DIFF_COUNTING_JASTROW_ORBITAL_H

#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/Jastrow/CountingRegion.h"

#include "QMCWaveFunctions/DiffWaveFunctionComponent.h"

namespace qmcplusplus{

template <class RegionType>
class DiffCountingJastrowOrbital: public DiffWaveFunctionComponent
{
  // variables

public:
  // constructor
  DiffCountingJastrowOrbital(ParticleSet& P)
  {
  }

  // extended checks and array resizing
  void initialize()
  {
  }
  
  void checkOutVariables(const opt_variables_type& optvars)
  {
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
    // 'prepare internal data for a new particle set'
  }

  void resetParameters(const opt_variables_type& optvars)
  {
  }

  void evaluateDerivRatios(ParticleSet& VP, const opt_variables_type& optvars, Matrix<ValueType>& dratios)
  {
  }

  void evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars, 
    std::vector<RealType>& dlogpsi, std::vector<RealType>& dhpsioverpsi)
  {
  }

  bool addRegion(RegionType* CR);

  bool addDebug(int seqlen, int period);

};

}
#endif
