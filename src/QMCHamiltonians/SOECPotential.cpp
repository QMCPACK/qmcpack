//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
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
    : myRNG(&Random), IonConfig(ions), Psi(psi), Peln(els), ElecNeighborIons(els), IonNeighborElecs(ions)
{
  //set_energy_domain(potential);
  //two_body_quantum_domain(ions,els);
  myTableIndex = els.addTable(ions, DT_SOA_PREFERRED);
  NumIons      = ions.getTotalNum();
  PP.resize(NumIons, nullptr);
  PPset.resize(IonConfig.getSpeciesSet().getTotalNum(), 0);
  //UpdateMode.set(NONLOCAL,1);
}

//destructor
SOECPotential::~SOECPotential() { delete_iter(PPset.begin(), PPset.end()); }

void SOECPotential::resetTargetParticleSet(ParticleSet& P) {}

SOECPotential::Return_t SOECPotential::evaluate(ParticleSet& P) { return 0.0; }

OperatorBase* SOECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi) { return 0; }

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
