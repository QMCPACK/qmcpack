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

#include "Particle/DistanceTable.h"
#include "SOECPotential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
/** constructor
 *\param ionic positions
 *\param els electronic poitions
 *\param psi Trial wave function
*/
SOECPotential::SOECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi)
    : myRNG_(nullptr), IonConfig_(ions), Psi_(psi), Peln_(els), ElecNeighborIons_(els), IonNeighborElecs_(ions)
{
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ions, els);
  myTableIndex_ = els.addTable(ions);
  NumIons_      = ions.getTotalNum();
  PP_.resize(NumIons_, nullptr);
  PPset_.resize(IonConfig_.getSpeciesSet().getTotalNum());
}

void SOECPotential::resetTargetParticleSet(ParticleSet& P) {}

SOECPotential::Return_t SOECPotential::evaluate(ParticleSet& P)
{
  value_ = 0.0;
  for (int ipp = 0; ipp < PPset_.size(); ipp++)
    if (PPset_[ipp])
      PPset_[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*myRNG_));

  const auto& myTable = P.getDistTableAB(myTableIndex_);
  for (int iat = 0; iat < NumIons_; iat++)
    IonNeighborElecs_.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons_.getNeighborList(jel).clear();

  for (int jel = 0; jel < P.getTotalNum(); jel++)
  {
    const auto& dist               = myTable.getDistRow(jel);
    const auto& displ              = myTable.getDisplRow(jel);
    std::vector<int>& NeighborIons = ElecNeighborIons_.getNeighborList(jel);
    for (int iat = 0; iat < NumIons_; iat++)
      if (PP_[iat] != nullptr && dist[iat] < PP_[iat]->getRmax())
      {
        RealType pairpot = PP_[iat]->evaluateOne(P, iat, Psi_, jel, dist[iat], -displ[iat]);
        value_ += pairpot;
        NeighborIons.push_back(iat);
        IonNeighborElecs_.getNeighborList(iat).push_back(jel);
      }
  }
  return value_;
}

SOECPotential::Return_t SOECPotential::evaluateValueAndDerivatives(ParticleSet& P,
                                                                   const opt_variables_type& optvars,
                                                                   const Vector<ValueType>& dlogpsi,
                                                                   Vector<ValueType>& dhpsioverpsi)
{
  value_ = 0.0;
  for (int ipp = 0; ipp < PPset_.size(); ipp++)
    if (PPset_[ipp])
      PPset_[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*myRNG_));

  const auto& myTable = P.getDistTableAB(myTableIndex_);
  for (int jel = 0; jel < P.getTotalNum(); jel++)
  {
    const auto& dist  = myTable.getDistRow(jel);
    const auto& displ = myTable.getDisplRow(jel);
    for (int iat = 0; iat < NumIons_; iat++)
      if (PP_[iat] != nullptr && dist[iat] < PP_[iat]->getRmax())
        value_ += PP_[iat]->evaluateValueAndDerivatives(P, iat, Psi_, jel, dist[iat], -displ[iat], optvars, dlogpsi,
                                                        dhpsioverpsi);
  }
  return value_;
}

std::unique_ptr<OperatorBase> SOECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<SOECPotential> myclone = std::make_unique<SOECPotential>(IonConfig_, qp, psi);
  for (int ig = 0; ig < PPset_.size(); ++ig)
    if (PPset_[ig])
      myclone->addComponent(ig, std::unique_ptr<SOECPComponent>(PPset_[ig]->makeClone(qp)));
  return myclone;
}

void SOECPotential::addComponent(int groupID, std::unique_ptr<SOECPComponent>&& ppot)
{
  for (int iat = 0; iat < PP_.size(); iat++)
    if (IonConfig_.GroupID[iat] == groupID)
      PP_[iat] = ppot.get();
  PPset_[groupID] = std::move(ppot);
}

} // namespace qmcplusplus
