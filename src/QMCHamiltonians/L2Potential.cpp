//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "Particle/ParticleSet.h"
#include "DistanceTable.h"
#include "L2Potential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
L2Potential::L2Potential(const ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi) : IonConfig(ions)
{
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ions, els);
  NumIons      = ions.getTotalNum();
  myTableIndex = els.addTable(ions);
  size_t ns    = ions.getSpeciesSet().getTotalNum();
  PPset.resize(ns);
  PP.resize(NumIons, nullptr);
  psi_ref = &psi;
}

void L2Potential::resetTargetParticleSet(ParticleSet& P)
{
  int tid = P.addTable(IonConfig);
  if (tid != myTableIndex)
  {
    APP_ABORT("  L2Potential::resetTargetParticleSet found a different distance table index.");
  }
}


void L2Potential::add(int groupID, std::unique_ptr<L2RadialPotential>&& ppot)
{
  for (int iat = 0; iat < PP.size(); iat++)
    if (IonConfig.GroupID[iat] == groupID)
      PP[iat] = ppot.get();
  PPset[groupID] = std::move(ppot);
}


L2Potential::Return_t L2Potential::evaluate(ParticleSet& P)
{
  // compute the Hessian
  TrialWaveFunction::HessVector D2;
  // evaluateHessian gives the Hessian(log(Psi))
  psi_ref->evaluateHessian(P, D2);
  // add gradient terms to get (Hessian(Psi))/Psi instead
  const size_t N = P.getTotalNum();
  for (size_t n = 0; n < N; n++)
    for (size_t i = 0; i < DIM; i++)
      for (size_t j = 0; j < DIM; j++)
        D2[n](i, j) += P.G[n][i] * P.G[n][j];

  // compute v_L2(r)*L^2 for all electron-ion pairs
  const auto& d_table(P.getDistTableAB(myTableIndex));
  value_             = 0.0;
  const size_t Nelec = P.getTotalNum();
  for (size_t iel = 0; iel < Nelec; ++iel)
  {
    const auto Le    = P.L[iel];
    const auto& ge   = P.G[iel];
    const auto& D2e  = D2[iel];
    const auto& dist = d_table.getDistRow(iel);
    const auto& disp = d_table.getDisplRow(iel);
    Return_t esum    = 0.0;
    for (size_t iat = 0; iat < NumIons; ++iat)
    {
      RealType r = dist[iat];
      if (PP[iat] != nullptr && r < PP[iat]->rcut)
      {
        PosType rv = disp[iat]; //SoA rv is r_I-r_e
        RealType v = -r * r * std::real(Le);
        for (int i = 0; i < DIM; ++i)
          v += -r * r * std::real(ge[i] * ge[i]) - 2 * rv[i] * std::real(ge[i]);
        for (int i = 0; i < DIM; ++i)
          for (int j = 0; j < DIM; ++j)
            v += rv[i] * std::real(D2e(i, j)) * rv[j];
        esum += v * PP[iat]->evaluate(dist[iat]);
      }
    }
    value_ += esum;
  }
  return value_;
}


void L2Potential::evaluateDK(ParticleSet& P, int iel, TensorType& D, PosType& K)
{
  K = 0.0;
  D = 0.0;
  D.diagonal(1.0);

  const auto& d_table(P.getDistTableAB(myTableIndex));

  for (int iat = 0; iat < NumIons; iat++)
  {
    L2RadialPotential* ppot = PP[iat];
    if (ppot == nullptr)
      continue;
    RealType r = d_table.getTempDists()[iat];
    if (r < ppot->rcut)
    {
      PosType rv   = -1 * d_table.getTempDispls()[iat];
      RealType vL2 = ppot->evaluate(r);
      K += 2 * rv * vL2;
      for (int i = 0; i < DIM; ++i)
        D(i, i) += 2 * vL2 * r * r;
      for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
          D(i, j) -= 2 * vL2 * rv[i] * rv[j];
    }
  }
}


void L2Potential::evaluateD(ParticleSet& P, int iel, TensorType& D)
{
  D = 0.0;
  D.diagonal(1.0);

  const auto& d_table(P.getDistTableAB(myTableIndex));

  for (int iat = 0; iat < NumIons; iat++)
  {
    L2RadialPotential* ppot = PP[iat];
    if (ppot == nullptr)
      continue;
    RealType r = d_table.getTempDists()[iat];
    if (r < ppot->rcut)
    {
      PosType rv   = d_table.getTempDispls()[iat];
      RealType vL2 = ppot->evaluate(r);
      for (int i = 0; i < DIM; ++i)
        D(i, i) += 2 * vL2 * r * r;
      for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
          D(i, j) -= 2 * vL2 * rv[i] * rv[j];
    }
  }
}


std::unique_ptr<OperatorBase> L2Potential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<L2Potential> myclone = std::make_unique<L2Potential>(IonConfig, qp, psi);
  for (int ig = 0; ig < PPset.size(); ++ig)
  {
    if (PPset[ig])
    {
      myclone->add(ig, std::unique_ptr<L2RadialPotential>(PPset[ig]->makeClone()));
    }
  }
  return myclone;
}

} // namespace qmcplusplus
