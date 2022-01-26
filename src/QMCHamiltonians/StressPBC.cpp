//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Yubo Paul Yang, yyang173@illinois.edu, University of Illinois Urbana-Champaign
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "StressPBC.h"
#include "DistanceTable.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{
StressPBC::StressPBC(ParticleSet& ions, ParticleSet& elns, TrialWaveFunction& Psi0)
    : ForceBase(ions, elns),
      Psi(Psi0),
      PtclTarg(elns),
      PtclA(ions),
      ei_table_index(elns.addTable(ions)),
      ee_table_index(elns.addTable(elns)),
      ii_table_index(ions.addTable(ions)),
      firstTimeStress(true)
{
  ReportEngine PRE("StressPBC", "StressPBC");
  name_  = "StressPBC";
  prefix = "StressPBC";
  //This sets up the long range breakups.
  initBreakup(PtclTarg);
  stress_eI_const = 0.0;
  stress_ee_const = 0.0;
  if (firstTimeStress)
  { // calculate constants
    CalculateIonIonStress();
    firstTimeStress = false;
  }
  RealType vinv = -1. / PtclTarg.getLattice().Volume;
  app_log() << "\n====ion-ion stress ====\n" << stress_IonIon * vinv << std::endl;
  app_log() << "\n e-e const = " << stress_ee_const * vinv << std::endl;
  app_log() << "\n e-I const = " << stress_eI_const * vinv << std::endl;
}

void StressPBC::initBreakup(ParticleSet& P)
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
  NofSpeciesA.resize(NumSpeciesA);
  NofSpeciesB.resize(NumSpeciesB);
  for (int spec = 0; spec < NumSpeciesA; spec++)
  {
    Zspec[spec]       = tspeciesA(ChargeAttribIndxA, spec);
    NofSpeciesA[spec] = static_cast<int>(tspeciesA(MemberAttribIndxA, spec));
  }
  for (int spec = 0; spec < NumSpeciesB; spec++)
  {
    Qspec[spec]       = tspeciesB(ChargeAttribIndxB, spec);
    NofSpeciesB[spec] = static_cast<int>(tspeciesB(MemberAttribIndxB, spec));
  }

  for (int spec = 0; spec < NumSpeciesA; spec++) {}
  for (int iat = 0; iat < NptclA; iat++)
    Zat[iat] = Zspec[PtclA.GroupID[iat]];
  for (int iat = 0; iat < NptclB; iat++)
    Qat[iat] = Qspec[P.GroupID[iat]];

  AA = LRCoulombSingleton::getDerivHandler(P);
}

SymTensor<StressPBC::RealType, OHMMS_DIM> StressPBC::evaluateLR_AB(ParticleSet& P)
{
  SymTensor<RealType, OHMMS_DIM> res = 0.0;
  const StructFact& RhoKA(PtclA.getSK());
  const StructFact& RhoKB(P.getSK());

  for (int i = 0; i < NumSpeciesA; i++)
  {
    SymTensor<RealType, OHMMS_DIM> esum;
    esum = 0.0;
    for (int j = 0; j < NumSpeciesB; j++)
    {
#if defined(USE_REAL_STRUCT_FACTOR)
      esum += Qspec[j] *
          AA->evaluateStress(P.getSimulationCell().getKLists().kshell, RhoKA.rhok_r[i], RhoKA.rhok_i[i], RhoKB.rhok_r[j],
                             RhoKB.rhok_i[j]);
#else
      esum += Qspec[j] * AA->evaluateStress(P.getSimulationCell().getKLists().kshell, RhoKA.rhok[i], RhoKB.rhok[j]);

#endif
    }
    res += Zspec[i] * esum;
  }

  return res;
}

SymTensor<StressPBC::RealType, OHMMS_DIM> StressPBC::evaluateSR_AB(ParticleSet& P)
{
  const auto& d_ab                   = P.getDistTableAB(ei_table_index);
  SymTensor<RealType, OHMMS_DIM> res = 0.0;
  //Loop over distinct eln-ion pairs
  for (int jpart = 0; jpart < NptclB; jpart++)
  {
    const auto& drijs = d_ab.getDisplRow(jpart);
    const auto& rijs  = d_ab.getDistRow(jpart);
    const RealType e  = Qat[jpart];
    for (int iat = 0; iat < NptclA; iat++)
    {
      res += Zat[iat] * e * AA->evaluateSR_dstrain(drijs[iat], rijs[iat]);
    }
  }
  return res;
}

SymTensor<StressPBC::RealType, OHMMS_DIM> StressPBC::evaluateSR_AA(ParticleSet& P, int itabSelf)
{
  const auto& d_aa = P.getDistTableAA(itabSelf);

  SymTensor<RealType, OHMMS_DIM> stress_aa;
  for (int ipart = 0; ipart < NptclB; ipart++)
  {
    SymTensor<RealType, OHMMS_DIM> esum = 0.0;
    const auto& drijs                   = d_aa.getDisplRow(ipart);
    const auto& rijs                    = d_aa.getDistRow(ipart);
    for (int jpart = 0; jpart < ipart; jpart++)
    {
      esum += P.Z[jpart] * AA->evaluateSR_dstrain(drijs[jpart], rijs[jpart]);
    }
    stress_aa += P.Z[ipart] * esum;
  }

  return stress_aa;
}

SymTensor<StressPBC::RealType, OHMMS_DIM> StressPBC::evaluateLR_AA(ParticleSet& P)
{
  int NumSpecies = P.getSpeciesSet().TotalNum;
  SymTensor<RealType, OHMMS_DIM> stress_aa;
  const StructFact& PtclRhoK(P.getSK());
  int ChargeAttribIndx = P.getSpeciesSet().getAttribute("charge");
  int MemberAttribIndx = P.getSpeciesSet().getAttribute("membersize");

  std::vector<int> NofSpecies;
  std::vector<int> Zmyspec;
  NofSpecies.resize(NumSpecies);
  Zmyspec.resize(NumSpecies);

  for (int spec = 0; spec < NumSpecies; spec++)
  {
    Zmyspec[spec]    = P.getSpeciesSet()(ChargeAttribIndx, spec);
    NofSpecies[spec] = static_cast<int>(P.getSpeciesSet()(MemberAttribIndx, spec));
  }

  SymTensor<RealType, OHMMS_DIM> temp;
  for (int spec1 = 0; spec1 < NumSpecies; spec1++)
  {
    RealType Z1 = Zmyspec[spec1];
    for (int spec2 = spec1; spec2 < NumSpecies; spec2++)
    {
#if !defined(USE_REAL_STRUCT_FACTOR)
      SymTensor<RealType, OHMMS_DIM> temp =
          AA->evaluateStress(P.getSimulationCell().getKLists().kshell, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
#else
      SymTensor<RealType, OHMMS_DIM> temp =
          AA->evaluateStress(P.getSimulationCell().getKLists().kshell, PtclRhoK.rhok_r[spec1], PtclRhoK.rhok_i[spec1],
                             PtclRhoK.rhok_r[spec2], PtclRhoK.rhok_i[spec2]);
#endif
      if (spec2 == spec1)
        temp *= 0.5;
      stress_aa += Z1 * Zmyspec[spec2] * temp;
    } //spec2
  }   //spec1

  return stress_aa;
}

SymTensor<StressPBC::RealType, OHMMS_DIM> StressPBC::evalConsts_AB()
{
  int nelns = PtclTarg.getTotalNum();
  int nions = PtclA.getTotalNum();

  using mRealType = LRHandlerType::mRealType;

  SymTensor<mRealType, OHMMS_DIM> Consts = 0.0;
  SymTensor<mRealType, OHMMS_DIM> vs_k0  = AA->evaluateSR_k0_dstrain();
  mRealType v1; //single particle energy
  for (int i = 0; i < nelns; ++i)
  {
    v1 = 0.0;
    for (int s = 0; s < NumSpeciesA; s++)
      v1 += NofSpeciesA[s] * Zspec[s];
    Consts += -.5 * Qat[i] * vs_k0 * v1;
  }
  for (int i = 0; i < nions; ++i)
  {
    v1 = 0.0;
    for (int s = 0; s < NumSpeciesB; s++)
      v1 += NofSpeciesB[s] * Qspec[s];
    Consts += -.5 * Zat[i] * vs_k0 * v1;
  }

  return SymTensor<RealType, OHMMS_DIM>(Consts);
}

SymTensor<StressPBC::RealType, OHMMS_DIM> StressPBC::evalConsts_AA(ParticleSet& P)
{
  SymTensor<RealType, OHMMS_DIM> tmpconsts = 0.0; // constant term

  int NumSpecies = P.getSpeciesSet().TotalNum;
  int NumCenters = P.getTotalNum();
  RealType v1; //single particle energy

  int ChargeAttribIndx = P.getSpeciesSet().getAttribute("charge");
  int MemberAttribIndx = P.getSpeciesSet().getAttribute("membersize");

  std::vector<int> NofSpecies;
  std::vector<int> Zmyspec;
  NofSpecies.resize(NumSpecies);
  Zmyspec.resize(NumSpecies);

  for (int spec = 0; spec < NumSpecies; spec++)
  {
    Zmyspec[spec]    = P.getSpeciesSet()(ChargeAttribIndx, spec);
    NofSpecies[spec] = static_cast<int>(P.getSpeciesSet()(MemberAttribIndx, spec));
  }

  SymTensor<RealType, OHMMS_DIM> vl_r0 = AA->evaluateLR_r0_dstrain();
  for (int ipart = 0; ipart < NumCenters; ipart++)
  {
    tmpconsts += -.5 * P.Z[ipart] * P.Z[ipart] * vl_r0;
  }
  // if(report)
  app_log() << "   PBCAA self-interaction term \n" << tmpconsts << std::endl;
  //Compute Madelung constant: this is not correct for general cases
  SymTensor<RealType, OHMMS_DIM> MC0 = 0;
  for (int ks = 0; ks < AA->Fk.size(); ks++)
    MC0 += AA->dFk_dstrain[ks];
  MC0 = 0.5 * (MC0 - vl_r0);
  //Neutraling background term
  SymTensor<RealType, OHMMS_DIM> vs_k0 = AA->evaluateSR_k0_dstrain(); //v_s(k=0)
  for (int ipart = 0; ipart < NumCenters; ipart++)
  {
    v1 = 0.0;
    for (int spec = 0; spec < NumSpecies; spec++)
      v1 += -.5 * P.Z[ipart] * NofSpecies[spec] * Zmyspec[spec];
    tmpconsts += v1 * vs_k0;
  }
  // if(report)
  app_log() << "   PBCAA total constant \n" << tmpconsts << std::endl;
  return tmpconsts;
}

StressPBC::Return_t StressPBC::evaluate(ParticleSet& P)
{
  const RealType vinv(-1.0 / P.getLattice().Volume);
  stress     = 0.0;
  stress_ee  = 0.0;
  stress_ei  = 0.0;
  stress_kin = 0.0;

  stress_ei += vinv * evaluateLR_AB(PtclTarg);
  stress_ei += vinv * evaluateSR_AB(PtclTarg);
  stress_ei += vinv * stress_eI_const;

  stress_ee += vinv * evaluateLR_AA(PtclTarg);
  stress_ee += vinv * evaluateSR_AA(PtclTarg, ee_table_index);
  stress_ee += vinv * stress_ee_const;


  stress_kin += vinv * evaluateKineticSymTensor(P);

  stress = stress_ee + stress_ei + stress_kin;
  if (addionion)
    stress += vinv * stress_IonIon;

  return 0.0;
}

SymTensor<StressPBC::RealType, OHMMS_DIM> StressPBC::evaluateKineticSymTensor(ParticleSet& P)
{
  WaveFunctionComponent::HessVector grad_grad_psi;
  Psi.evaluateHessian(P, grad_grad_psi);
  SymTensor<RealType, OHMMS_DIM> kinetic_tensor;
  Tensor<ComplexType, OHMMS_DIM> complex_ktensor;

  for (int iat = 0; iat < P.getTotalNum(); iat++)
  {
    const RealType minv(1.0 / P.Mass[iat]);
    complex_ktensor += outerProduct(P.G[iat], P.G[iat]) * static_cast<ParticleSet::SingleParticleValue>(minv);
    complex_ktensor += grad_grad_psi[iat] * minv;
  }

  for (int i = 0; i < OHMMS_DIM; i++)
    for (int j = i; j < OHMMS_DIM; j++)
    {
      kinetic_tensor(i, j) = complex_ktensor(i, j).real();
    }
  return kinetic_tensor;
}

bool StressPBC::put(xmlNodePtr cur)
{
  std::string ionionforce("yes");
  OhmmsAttributeSet attr;
  attr.add(prefix, "name");
  attr.add(ionionforce, "addionion");
  attr.put(cur);
  addionion = (ionionforce == "yes") || (ionionforce == "true");
  app_log() << "add ion-ion stress = " << addionion << std::endl;
  return true;
}

std::unique_ptr<OperatorBase> StressPBC::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<StressPBC> tmp = std::make_unique<StressPBC>(PtclA, qp, psi);
  tmp->firstTimeStress           = firstTimeStress;
  tmp->stress_IonIon             = stress_IonIon;
  tmp->stress_ee_const           = stress_ee_const;
  tmp->stress_eI_const           = stress_eI_const;
  tmp->addionion                 = addionion;
  return tmp;
}
} // namespace qmcplusplus
