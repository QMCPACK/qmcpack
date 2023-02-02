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
template<typename T>
struct NLPPJob;

class SOECPotential : public OperatorBase
{
public:
  SOECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi);

  bool dependsOnWaveFunction() const override { return true; }
  std::string getClassName() const override { return "SOECPotential"; }
  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  Return_t evaluateValueAndDerivatives(ParticleSet& P,
                                       const opt_variables_type& optvars,
                                       const Vector<ValueType>& dlogpsi,
                                       Vector<ValueType>& dhpsioverpsi) override;

  void mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                   const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                   const RefVectorWithLeader<ParticleSet>& p_list) const override;

  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "SOECPotential: " << IonConfig_.getName();
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  void addComponent(int groupID, std::unique_ptr<SOECPComponent>&& pp);

  void setRandomGenerator(RandomGenerator* rng) override { myRNG_ = rng; }

protected:
  RandomGenerator* myRNG_;
  std::vector<SOECPComponent*> PP_;
  std::vector<std::unique_ptr<SOECPComponent>> PPset_;
  ParticleSet& IonConfig_;
  TrialWaveFunction& Psi_;

private:
  ///number of ions
  int NumIons_;
  ///index of distance table for ion-el pair
  int myTableIndex_;
  ///reference to the electrons
  ParticleSet& Peln_;
  ///neighborlist of electrons
  NeighborLists ElecNeighborIons_;
  ///neighborlist of ions
  NeighborLists IonNeighborElecs_;
  //job list for evaluation
  std::vector<std::vector<NLPPJob<RealType>>> sopp_jobs_;
};
} // namespace qmcplusplus

#endif
