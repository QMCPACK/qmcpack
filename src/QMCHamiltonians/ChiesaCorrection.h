//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim and Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  void addEnergy(MCWalkerConfiguration &W, vector<RealType> &LocalEnergy);
#endif


  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
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
