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
#include "QMCHamiltonians/L2Potential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
L2Potential::L2Potential(const ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi) : IonConfig(ions)
{
  set_energy_domain(potential);
  two_body_quantum_domain(ions, els);
  NumIons      = ions.getTotalNum();
  myTableIndex = els.addTable(ions, DT_SOA_PREFERRED);
  size_t ns    = ions.getSpeciesSet().getTotalNum();
  PPset.resize(ns, 0);
  PP.resize(NumIons, nullptr);
  psi_ref = &psi;
}


///destructor
L2Potential::~L2Potential() { delete_iter(PPset.begin(), PPset.end()); }


void L2Potential::resetTargetParticleSet(ParticleSet& P)
{
  int tid = P.addTable(IonConfig, DT_SOA_PREFERRED);
  if (tid != myTableIndex)
  {
    APP_ABORT("  L2Potential::resetTargetParticleSet found a different distance table index.");
  }
}


void L2Potential::add(int groupID, L2RadialPotential* ppot)
{
  PPset[groupID] = ppot;
  for (int iat = 0; iat < PP.size(); iat++)
    if (IonConfig.GroupID[iat] == groupID)
      PP[iat] = ppot;
}


L2Potential::Return_t L2Potential::evaluate(ParticleSet& P)
{
  // compute the Hessian
  TrialWaveFunction::HessVector_t D2;
  // evaluateHessian gives the Hessian(log(Psi))
  psi_ref->evaluateHessian(P, D2);
  // add gradient terms to get (Hessian(Psi))/Psi instead
  const size_t N = P.getTotalNum();
  for (size_t n = 0; n < N; n++)
    for (size_t i = 0; i < DIM; i++)
      for (size_t j = 0; j < DIM; j++)
        D2[n](i, j) += P.G[n][i] * P.G[n][j];

  // compute v_L2(r)*L^2 for all electron-ion pairs
  const DistanceTableData& d_table(P.getDistTable(myTableIndex));
  Value = 0.0;
  if (d_table.DTType == DT_SOA)
  {
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
      Value += esum;
    }
  }
  else
  {
#ifndef ENABLE_SOA
    for (int iat = 0; iat < NumIons; iat++)
    {
      L2RadialPotential* ppot = PP[iat];
      if (ppot == nullptr)
        continue;
      Return_t esum = 0.0;
      for (int nn = d_table.M[iat], n = 0; nn < d_table.M[iat + 1]; ++nn, ++n)
      {
        RealType r = d_table.r(nn);
        if (r < ppot->rcut)
        {
          PosType rv      = d_table.dr(nn); //AoS rv is r_e-r_I
          const auto Le   = P.L[n];
          const auto& ge  = P.G[n];
          const auto& D2e = D2[n];
          RealType v      = -r * r * std::real(Le);
          for (int i = 0; i < DIM; ++i)
            v += -r * r * std::real(ge[i] * ge[i]) + 2 * rv[i] * std::real(ge[i]);
          for (int i = 0; i < DIM; ++i)
            for (int j = 0; j < DIM; ++j)
              v += rv[i] * std::real(D2e(i, j)) * rv[j];
          esum += v * ppot->evaluate(r);
        }
      }
      Value += esum;
    }
#endif
  }
  return Value;
}


OperatorBase* L2Potential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  L2Potential* myclone = new L2Potential(IonConfig, qp, psi);
  for (int ig = 0; ig < PPset.size(); ++ig)
    if (PPset[ig])
    {
      L2RadialPotential* ppot = PPset[ig]->makeClone();
      myclone->add(ig, ppot);
    }
  return myclone;
}

} // namespace qmcplusplus
