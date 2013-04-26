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

#ifndef QMCPLUSPLUS_PULAY_FOCE_H
#define QMCPLUSPLUS_PULAY_FOCE_H

#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/ForceBase.h"

namespace qmcplusplus
{

struct PulayForce : public QMCHamiltonianBase, public ForceBase
{
  ParticleSet& Ions;
  ParticleSet& Electrons;
  TrialWaveFunction& Psi;

  vector<RealType> WarpNorm;

  ParticleSet::ParticlePos_t GradLogPsi, EGradLogPsi;

  PulayForce(ParticleSet& ions, ParticleSet& elns,
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
    PulayForce *myClone = new PulayForce (Ions, qp, psi);
    myClone->FirstForceIndex = FirstForceIndex;
    return myClone;
  }


  inline RealType
  WarpFunction (RealType r)
  {
    return 1.0/(r*r*r*r);
  }

  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist);

  void setParticlePropertyList(PropertySetType& plist, int offset);

  void registerObservables(vector<observable_helper*>& h5list,
                           hid_t gid) const;
};

}
#endif
