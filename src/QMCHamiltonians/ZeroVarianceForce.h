//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Ken Esler
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#ifndef QMCPLUSPLUS_ZERO_VARIANCE_FORCE_H
#define QMCPLUSPLUS_ZERO_VARIANCE_FORCE_H

#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/ForceBase.h"

namespace qmcplusplus
{

struct ZeroVarianceForce : public QMCHamiltonianBase, public ForceBase
{
  ParticleSet& Ions;
  ParticleSet& Electrons;
  TrialWaveFunction& Psi;

  ParticleSet::ParticlePos_t F_ZV1, F_ZV2;
  TinyVector<ParticleSet::ParticleGradient_t,OHMMS_DIM>  grad_grad_psi;
  TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad_psi;

  ZeroVarianceForce(ParticleSet& ions, ParticleSet& elns,
                    TrialWaveFunction &psi);

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    ZeroVarianceForce *myClone = new ZeroVarianceForce (Ions, qp, psi);
    myClone->FirstForceIndex = FirstForceIndex;
    return myClone;
  }

  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist);

  void setParticlePropertyList(PropertySetType& plist, int offset);

  void registerObservables(vector<observable_helper*>& h5list,
                           hid_t gid) const;
};

}
#endif
