//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

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

  std::vector<RealType> WarpNorm;

  ParticleSet::ParticlePos_t GradLogPsi, EGradLogPsi;

  PulayForce(ParticleSet& ions, ParticleSet& elns,
             TrialWaveFunction &psi);

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
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

  void registerObservables(std::vector<observable_helper*>& h5list,
                           hid_t gid) const;
};

}
#endif
