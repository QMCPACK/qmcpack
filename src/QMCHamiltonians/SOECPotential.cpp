//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories 
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/SOECPotential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
/** constructor
 *\param ionic positions
 *\param els electronic poitions
 *\param psi Trial wave function
*/
SOECPotential::SOECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi)
    : myRNG(nullptr), IonConfig(ions), Psi(psi), Peln(els), ElecNeighborIons(els), IonNeighborElecs(ions)
{
  set_energy_domain(potential);
  two_body_quantum_domain(ions, els);
  myTableIndex = els.addTable(ions, DT_SOA_PREFERRED);
  NumIons      = ions.getTotalNum();
  PP.resize(NumIons, nullptr);
  PPset.resize(IonConfig.getSpeciesSet().getTotalNum(), 0);
}

//destructor
SOECPotential::~SOECPotential() { delete_iter(PPset.begin(), PPset.end()); }

void SOECPotential::resetTargetParticleSet(ParticleSet& P) {}

SOECPotential::Return_t SOECPotential::evaluate(ParticleSet& P)
{
  Value = 0.0;
  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      PPset[ipp]->randomize_grid(*myRNG);
  const auto& myTable = P.getDistTable(myTableIndex);
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  if (myTable.DTType == DT_SOA)
  {
    for (int jel = 0; jel < P.getTotalNum(); jel++)
    {
      const auto& dist               = myTable.getDistRow(jel);
      const auto& displ              = myTable.getDisplRow(jel);
      std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
      for (int iat = 0; iat < NumIons; iat++)
        if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
        {
          RealType pairpot = PP[iat]->evaluateOne(P, iat, Psi, jel, dist[iat], -displ[iat]);
          Value += pairpot;
          NeighborIons.push_back(iat);
          IonNeighborElecs.getNeighborList(iat).push_back(jel);
        }
    }
  }
  else
  {
    APP_ABORT("SOECPotential::evaluate(): AOS is deprecated. Distance tables must be SOA\n");
  }
  return Value;
}

OperatorBase* SOECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  SOECPotential* myclone = new SOECPotential(IonConfig, qp, psi);
  for (int ig = 0; ig < PPset.size(); ++ig)
  {
    if (PPset[ig])
    {
      SOECPComponent* ppot = PPset[ig]->makeClone(qp);
      myclone->addComponent(ig, ppot);
    }
  }
  return myclone;
}

void SOECPotential::addComponent(int groupID, SOECPComponent* ppot)
{
  for (int iat = 0; iat < PP.size(); iat++)
    if (IonConfig.GroupID[iat] == groupID)
    {
      PP[iat] = ppot;
    }
  PPset[groupID] = ppot;
}

} // namespace qmcplusplus
