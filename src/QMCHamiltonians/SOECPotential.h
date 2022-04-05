//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SO_ECPOTENTIAL_H
#define QMCPLUSPLUS_SO_ECPOTENTIAL_H

#include "QMCHamiltonians/SOECPComponent.h"
#include "Particle/NeighborLists.h"

namespace qmcplusplus
{
class SOECPotential : public OperatorBase
{
public:
  SOECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi);

  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "SOECPotential: " << IonConfig.getName();
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  void addComponent(int groupID, std::unique_ptr<SOECPComponent>&& pp);

  void setRandomGenerator(RandomGenerator* rng) override { myRNG = rng; }

protected:
  RandomGenerator* myRNG;
  std::vector<SOECPComponent*> PP;
  std::vector<std::unique_ptr<SOECPComponent>> PPset;
  ParticleSet& IonConfig;
  TrialWaveFunction& Psi;

private:
  ///number of ions
  int NumIons;
  ///index of distance table for ion-el pair
  int myTableIndex;
  ///reference to the electrons
  ParticleSet& Peln;
  ///neighborlist of electrons
  NeighborLists ElecNeighborIons;
  ///neighborlist of ions
  NeighborLists IonNeighborElecs;
};
} // namespace qmcplusplus

#endif
