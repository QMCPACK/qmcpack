//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "LocalECPotential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
LocalECPotential::LocalECPotential(const ParticleSet& ions, ParticleSet& els) : IonConfig(ions), Peln(els), Pion(ions)
{
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ions, els);
  NumIons      = ions.getTotalNum();
  myTableIndex = els.addTable(ions);
  //allocate null
  PPset.resize(ions.getSpeciesSet().getTotalNum());
  PP.resize(NumIons, nullptr);
  Zeff.resize(NumIons, 0.0);
  gZeff.resize(ions.getSpeciesSet().getTotalNum(), 0);
}

void LocalECPotential::resetTargetParticleSet(ParticleSet& P)
{
  int tid = P.addTable(IonConfig);
  if (tid != myTableIndex)
  {
    APP_ABORT("  LocalECPotential::resetTargetParticleSet found a different distance table index.");
  }
}

void LocalECPotential::add(int groupID, std::unique_ptr<RadialPotentialType>&& ppot, RealType z)
{
  for (int iat = 0; iat < PP.size(); iat++)
  {
    if (IonConfig.GroupID[iat] == groupID)
    {
      PP[iat]   = ppot.get();
      Zeff[iat] = z;
    }
  }
  PPset[groupID] = std::move(ppot);
  gZeff[groupID] = z;
}

#if !defined(REMOVE_TRACEMANAGER)
void LocalECPotential::contributeParticleQuantities() { request_.contribute_array(name_); }

void LocalECPotential::checkoutParticleQuantities(TraceManager& tm)
{
  streaming_particles_ = request_.streaming_array(name_);
  if (streaming_particles_)
  {
    Ve_sample = tm.checkout_real<1>(name_, Peln);
    Vi_sample = tm.checkout_real<1>(name_, Pion);
  }
}

void LocalECPotential::deleteParticleQuantities()
{
  if (streaming_particles_)
  {
    delete Ve_sample;
    delete Vi_sample;
  }
}
#endif


LocalECPotential::Return_t LocalECPotential::evaluate(ParticleSet& P)
{
#if !defined(REMOVE_TRACEMANAGER)
  if (streaming_particles_)
    value_ = evaluate_sp(P);
  else
#endif
  {
    const auto& d_table(P.getDistTableAB(myTableIndex));
    value_             = 0.0;
    const size_t Nelec = P.getTotalNum();
    for (size_t iel = 0; iel < Nelec; ++iel)
    {
      const auto& dist = d_table.getDistRow(iel);
      Return_t esum(0);
      for (size_t iat = 0; iat < NumIons; ++iat)
        if (PP[iat] != nullptr)
          esum += PP[iat]->splint(dist[iat]) * Zeff[iat] / dist[iat];
      value_ -= esum;
    }
  }
  return value_;
}

LocalECPotential::Return_t LocalECPotential::evaluateWithIonDerivs(ParticleSet& P,
                                                                   ParticleSet& ions,
                                                                   TrialWaveFunction& psi,
                                                                   ParticleSet::ParticlePos_t& hf_terms,
                                                                   ParticleSet::ParticlePos_t& pulay_terms)
{
  const auto& d_table(P.getDistTableAB(myTableIndex));
  value_             = 0.0;
  const size_t Nelec = P.getTotalNum();
  for (size_t iel = 0; iel < Nelec; ++iel)
  {
    const auto& dist = d_table.getDistRow(iel);
    const auto& dr   = d_table.getDisplRow(iel);
    Return_t esum(0);
    //value, radial derivative, and 2nd derivative of spline r*V.
    RealType v(0.0), dv(0.0), d2v(0.0);
    //radial derivative dV/dr
    RealType dvdr(0.0);
    RealType rinv(1.0);
    for (size_t iat = 0; iat < NumIons; ++iat)
    {
      if (PP[iat] != nullptr)
      {
        rinv = 1.0 / dist[iat];
        v    = PP[iat]->splint(dist[iat], dv, d2v);
        dvdr = -Zeff[iat] * (dv - v * rinv) * rinv; //the minus is because of charge of electron.
        hf_terms[iat] += dvdr * dr[iat] * rinv;
        esum += -Zeff[iat] * v * rinv;
      }
    }
    value_ += esum;
  }
  return value_;
}

#if !defined(REMOVE_TRACEMANAGER)
LocalECPotential::Return_t LocalECPotential::evaluate_sp(ParticleSet& P)
{
  const auto& d_table(P.getDistTableAB(myTableIndex));
  value_                      = 0.0;
  Array<RealType, 1>& Ve_samp = *Ve_sample;
  Array<RealType, 1>& Vi_samp = *Vi_sample;
  Ve_samp                     = 0.0;
  Vi_samp                     = 0.0;

  const size_t Nelec = P.getTotalNum();
  for (size_t iel = 0; iel < Nelec; ++iel)
  {
    const auto& dist = d_table.getDistRow(iel);
    Return_t esum(0), pairpot;
    for (size_t iat = 0; iat < NumIons; ++iat)
      if (PP[iat] != nullptr)
      {
        pairpot = -0.5 * PP[iat]->splint(dist[iat]) * Zeff[iat] / dist[iat];
        Vi_samp(iat) += pairpot;
        Ve_samp(iel) += pairpot;
        esum += pairpot;
      }
    value_ += esum;
  }
  value_ *= 2.0;

#if defined(TRACE_CHECK)
  RealType Vnow  = Value;
  RealType Visum = Vi_samp.sum();
  RealType Vesum = Ve_samp.sum();
  RealType Vsum  = Vesum + Visum;
  RealType Vorig = evaluate_orig(P);
  if (std::abs(Vsum - Vnow) > TraceManager::trace_tol)
  {
    app_log() << "accumtest: LocalECPotential::evaluate()" << std::endl;
    app_log() << "accumtest:   tot:" << Vnow << std::endl;
    app_log() << "accumtest:   sum:" << Vsum << std::endl;
    APP_ABORT("Trace check failed");
  }
  if (std::abs(Vesum - Visum) > TraceManager::trace_tol)
  {
    app_log() << "sharetest: LocalECPotential::evaluate()" << std::endl;
    app_log() << "sharetest:   e share:" << Vesum << std::endl;
    app_log() << "sharetest:   i share:" << Visum << std::endl;
    APP_ABORT("Trace check failed");
  }
  if (std::abs(Vorig - Vnow) > TraceManager::trace_tol)
  {
    app_log() << "versiontest: LocalECPotential::evaluate()" << std::endl;
    app_log() << "versiontest:   orig:" << Vorig << std::endl;
    app_log() << "versiontest:    mod:" << Vnow << std::endl;
    APP_ABORT("Trace check failed");
  }
#endif
  return value_;
}
#endif


LocalECPotential::Return_t LocalECPotential::evaluate_orig(ParticleSet& P)
{
  const auto& d_table(P.getDistTableAB(myTableIndex));
  value_             = 0.0;
  const size_t Nelec = P.getTotalNum();
  for (size_t iel = 0; iel < Nelec; ++iel)
  {
    const auto& dist = d_table.getDistRow(iel);
    Return_t esum(0);
    for (size_t iat = 0; iat < NumIons; ++iat)
      if (PP[iat] != nullptr)
        esum += PP[iat]->splint(dist[iat]) * Zeff[iat] / dist[iat];
    value_ -= esum;
  }
  return value_;
}

std::unique_ptr<OperatorBase> LocalECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<LocalECPotential> myclone = std::make_unique<LocalECPotential>(IonConfig, qp);

  for (int ig = 0; ig < PPset.size(); ++ig)
  {
    if (PPset[ig])
    {
      myclone->add(ig, std::unique_ptr<RadialPotentialType>(PPset[ig]->makeClone()), gZeff[ig]);
    }
  }
  return myclone;
}
} // namespace qmcplusplus
