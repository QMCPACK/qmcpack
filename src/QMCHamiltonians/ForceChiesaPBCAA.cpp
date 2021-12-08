//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "ForceChiesaPBCAA.h"
#include "Particle/DistanceTable.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
ForceChiesaPBCAA::ForceChiesaPBCAA(ParticleSet& ions, ParticleSet& elns, bool firsttime)
    : ForceBase(ions, elns),
      PtclA(ions),
      first_time(firsttime),
      d_aa_ID(ions.addTable(ions)),
      d_ei_ID(elns.addTable(ions))
{
  ReportEngine PRE("ForceChiesaPBCAA", "ForceChiesaPBCAA");
  name_  = "Chiesa_Force_Base_PBCAB";
  prefix = "FChiesaPBC";
  //Defaults for the chiesa S-wave polynomial filtering.
  Rcut          = 0.4;
  m_exp         = 2;
  N_basis       = 4;
  forces        = 0.0;
  forces_IonIon = 0.0;
  ions.turnOnPerParticleSK();
  //This sets up the long range breakups.
  initBreakup(elns);
  if (first_time == true)
  { // calculate ion-ion forces and print in output
    ions.update();
    evaluateLR_AA();
    evaluateSR_AA();
    app_log() << "IonIon Force" << std::endl;
    app_log() << forces_IonIon << std::endl;
    first_time = false;
  }
}

void ForceChiesaPBCAA::InitMatrix()
{
  Sinv.resize(N_basis, N_basis);
  h.resize(N_basis);
  c.resize(N_basis);
  for (int k = 0; k < N_basis; k++)
  {
    h[k] = std::pow(Rcut, (k + 2)) / static_cast<RealType>(k + 2);
    for (int j = 0; j < N_basis; j++)
    {
      Sinv(k, j) = std::pow(Rcut, (m_exp + k + j + 3)) / static_cast<RealType>(m_exp + k + j + 3);
    }
  }
  // in Numerics/DeterminantOperators.h
  invert_matrix(Sinv, false);
  // in Numerics/MatrixOperators.h
  MatrixOperators::product(Sinv, h.data(), c.data());
}

void ForceChiesaPBCAA::initBreakup(ParticleSet& P)
{
  SpeciesSet& tspeciesA(PtclA.getSpeciesSet());
  SpeciesSet& tspeciesB(P.getSpeciesSet());
  int ChargeAttribIndxA = tspeciesA.addAttribute("charge");
  int MemberAttribIndxA = tspeciesA.addAttribute("membersize");
  int ChargeAttribIndxB = tspeciesB.addAttribute("charge");
  int MemberAttribIndxB = tspeciesB.addAttribute("membersize");
  NptclA                = PtclA.getTotalNum();
  NptclB                = P.getTotalNum();
  NumSpeciesA           = tspeciesA.TotalNum;
  NumSpeciesB           = tspeciesB.TotalNum;
  //Store information about charges and number of each species
  Zat.resize(NptclA);
  Zspec.resize(NumSpeciesA);
  Qat.resize(NptclB);
  Qspec.resize(NumSpeciesB);
  for (int spec = 0; spec < NumSpeciesA; spec++)
  {
    Zspec[spec] = tspeciesA(ChargeAttribIndxA, spec);
  }
  for (int spec = 0; spec < NumSpeciesB; spec++)
  {
    Qspec[spec] = tspeciesB(ChargeAttribIndxB, spec);
  }
  for (int iat = 0; iat < NptclA; iat++)
    Zat[iat] = Zspec[PtclA.GroupID[iat]];
  for (int iat = 0; iat < NptclB; iat++)
    Qat[iat] = Qspec[P.GroupID[iat]];
  dAB = LRCoulombSingleton::getDerivHandler(P);
}

void ForceChiesaPBCAA::evaluateLR(ParticleSet& P)
{
  std::vector<TinyVector<RealType, DIM>> grad(PtclA.getTotalNum());
  for (int j = 0; j < NumSpeciesB; j++)
  {
    for (int iat = 0; iat < grad.size(); iat++)
      grad[iat] = TinyVector<RealType, DIM>(0.0, 0.0, 0.0);
    dAB->evaluateGrad(PtclA, P, j, Zat, grad);
    for (int iat = 0; iat < grad.size(); iat++)
    {
      forces[iat] += Qspec[j] * grad[iat];
    }
  } // electron species
}

void ForceChiesaPBCAA::evaluateSR(ParticleSet& P)
{
  const auto& d_ab(P.getDistTableAB(d_ei_ID));
  for (size_t jat = 0; jat < NptclB; ++jat)
  {
    const auto& dist  = d_ab.getDistRow(jat);
    const auto& displ = d_ab.getDisplRow(jat);
    for (size_t iat = 0; iat < NptclA; ++iat)
    {
      const RealType r    = dist[iat];
      const RealType rinv = RealType(1) / r;
      RealType g_f        = g_filter(r);
      RealType V          = -dAB->srDf(r, rinv);
      PosType drhat       = rinv * displ[iat];
      forces[iat] += g_f * Zat[iat] * Qat[jat] * V * drhat;
    }
  }
}

void ForceChiesaPBCAA::evaluateSR_AA()
{
  const auto& d_aa(PtclA.getDistTableAA(d_aa_ID));
  for (size_t ipart = 1; ipart < NptclA; ipart++)
  {
    const auto& dist  = d_aa.getDistRow(ipart);
    const auto& displ = d_aa.getDisplRow(ipart);
    for (size_t jpart = 0; jpart < ipart; ++jpart)
    {
      RealType V   = -dAB->srDf(dist[jpart], RealType(1) / dist[jpart]);
      PosType grad = -Zat[jpart] * Zat[ipart] * V / dist[jpart] * displ[jpart];
      forces_IonIon[ipart] += grad;
      forces_IonIon[jpart] -= grad;
    }
  }
}

void ForceChiesaPBCAA::evaluateLR_AA()
{
  std::vector<TinyVector<RealType, DIM>> grad(PtclA.getTotalNum());
  for (int spec2 = 0; spec2 < NumSpeciesA; spec2++)
  {
    RealType Z2 = Zspec[spec2];
    for (int iat = 0; iat < grad.size(); iat++)
      grad[iat] = TinyVector<RealType, DIM>(0.0);
    dAB->evaluateGrad(PtclA, PtclA, spec2, Zat, grad);

    for (int iat = 0; iat < grad.size(); iat++)
    {
      forces_IonIon[iat] += Z2 * grad[iat];
    }
  } //spec2
}


ForceChiesaPBCAA::Return_t ForceChiesaPBCAA::evaluate(ParticleSet& P)
{
  forces = 0.0;
  evaluateLR(P);
  evaluateSR(P);
  if (addionion == true)
    forces = forces + forces_IonIon;
  return 0.0;
}

ForceChiesaPBCAA::Return_t ForceChiesaPBCAA::g_filter(RealType r)
{
  if (r >= Rcut)
  {
    return 1.0;
  }
  else
  {
    RealType g_q = 0.0;
    for (int q = 0; q < N_basis; q++)
    {
      g_q += c[q] * std::pow(r, m_exp + q + 1);
    }

    return g_q;
  }
}

bool ForceChiesaPBCAA::put(xmlNodePtr cur)
{
  std::string ionionforce("yes");
  OhmmsAttributeSet attr;
  attr.add(prefix, "name");
  attr.add(ionionforce, "addionion");
  attr.put(cur);
  addionion = (ionionforce == "yes") || (ionionforce == "true");
  app_log() << "ionionforce = " << ionionforce << std::endl;
  app_log() << "addionion=" << addionion << std::endl;
  ParameterSet fcep_param_set;
  fcep_param_set.add(Rcut, "rcut");
  fcep_param_set.add(N_basis, "nbasis");
  fcep_param_set.add(m_exp, "weight_exp");
  fcep_param_set.put(cur);
  app_log() << "    ForceChiesaPBCAA Parameters" << std::endl;
  app_log() << "        ForceChiesaPBCAA::Rcut=" << Rcut << std::endl;
  app_log() << "        ForceChiesaPBCAA::N_basis=" << N_basis << std::endl;
  app_log() << "        ForceChiesaPBCAA::m_exp=" << m_exp << std::endl;
  InitMatrix();
  return true;
}

void ForceChiesaPBCAA::resetTargetParticleSet(ParticleSet& P) { dAB->resetTargetParticleSet(P); }

void ForceChiesaPBCAA::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = plist.add(name_.c_str());
  addObservablesF(plist);
}

std::unique_ptr<OperatorBase> ForceChiesaPBCAA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<ForceChiesaPBCAA> tmp = std::make_unique<ForceChiesaPBCAA>(PtclA, qp, false);
  tmp->Rcut                             = Rcut;    // parameter: radial distance within which estimator is used
  tmp->m_exp                            = m_exp;   // parameter: exponent in polynomial fit
  tmp->N_basis                          = N_basis; // parameter: size of polynomial basis set
  tmp->Sinv.resize(N_basis, N_basis);
  tmp->Sinv = Sinv; // terms in fitting polynomial
  tmp->h.resize(N_basis);
  tmp->h = h; // terms in fitting polynomial
  tmp->c.resize(N_basis);
  tmp->c             = c; // polynomial coefficients
  tmp->addionion     = addionion;
  tmp->forces_IonIon = forces_IonIon;
  tmp->initBreakup(qp);
  return tmp;
}
} // namespace qmcplusplus
