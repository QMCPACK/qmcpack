//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_CHIESA_CORRECTION_H
#define QMCPLUSPLUS_CHIESA_CORRECTION_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{

class ChiesaCorrection : public QMCHamiltonianBase
{
private:
  const TrialWaveFunction &psi_ref;
  ParticleSet &ptcl_ref;

public:
  ChiesaCorrection (ParticleSet& ptcl, const TrialWaveFunction &psi) :
    psi_ref(psi), ptcl_ref(ptcl)
  {
  }

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

#ifdef QMC_CUDA
  void addEnergy(MCWalkerConfiguration &W, std::vector<RealType> &LocalEnergy);
#endif


  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  bool put(xmlNodePtr cur);

  bool get(std::ostream& os) const
  {
    os << "Chiesa correction: " << ptcl_ref.getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
};

}

#endif
