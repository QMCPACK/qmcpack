//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file HamiltonianFactory.cpp
 *@brief Manage instantiation of Coulomb-related potentials
 */
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/CoulombPBCAA.h"
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "QMCHamiltonians/ForceChiesaPBCAA.h"
#if OHMMS_DIM == 3
#include "QMCHamiltonians/ECPotentialBuilder.h"
#include "QMCHamiltonians/ForceBase.h"
#include "QMCHamiltonians/ForceCeperley.h"

#include "QMCHamiltonians/ACForce.h"
#if defined(HAVE_LIBFFTW)
#include "QMCHamiltonians/MPC.h"
#endif
#endif
#include "OhmmsData/AttributeSet.h"

#ifdef QMC_CUDA
#include "QMCHamiltonians/CoulombPBCAA_CUDA.h"
#include "QMCHamiltonians/CoulombPBCAB_CUDA.h"
#include "QMCHamiltonians/CoulombPotential_CUDA.h"
#include "QMCHamiltonians/MPC_CUDA.h"
#endif

//#include <iostream>
namespace qmcplusplus
{
void HamiltonianFactory::addMPCPotential(xmlNodePtr cur, bool isphysical)
{
#if OHMMS_DIM == 3 && defined(HAVE_LIBFFTW)
  std::string a("e"), title("MPC"), physical("no");
  OhmmsAttributeSet hAttrib;
  double cutoff = 30.0;
  hAttrib.add(title, "id");
  hAttrib.add(title, "name");
  hAttrib.add(cutoff, "cutoff");
  hAttrib.add(physical, "physical");
  hAttrib.put(cur);
  renameProperty(a);
  isphysical = (physical == "yes" || physical == "true");
#ifdef QMC_CUDA
  std::unique_ptr<MPC_CUDA> mpc = std::make_unique<MPC_CUDA>(targetPtcl, cutoff);
#else
  std::unique_ptr<MPC> mpc = std::make_unique<MPC>(targetPtcl, cutoff);
#endif
  targetH->addOperator(std::move(mpc), "MPC", isphysical);
#else
  APP_ABORT(
      "HamiltonianFactory::addMPCPotential MPC is disabled because FFTW3 was not found during the build process.");
#endif // defined(HAVE_LIBFFTW)
}


void HamiltonianFactory::addCoulombPotential(xmlNodePtr cur)
{
  using Return_t = QMCHamiltonian::FullPrecRealType;
  std::string targetInp(targetPtcl.getName());
  std::string sourceInp(targetPtcl.getName());
  std::string title("ElecElec"), pbc("yes");
  std::string forces("no");
  std::string use_gpu;
  bool physical = true;
  OhmmsAttributeSet hAttrib;
  hAttrib.add(title, "id");
  hAttrib.add(title, "name");
  hAttrib.add(targetInp, "target");
  hAttrib.add(sourceInp, "source");
  hAttrib.add(pbc, "pbc");
  hAttrib.add(physical, "physical");
  hAttrib.add(forces, "forces");
  hAttrib.add(use_gpu, "gpu", {"", "yes", "no"});
  hAttrib.put(cur);
  const bool applyPBC = (PBCType && pbc == "yes");
  const bool doForces = (forces == "yes") || (forces == "true");
  ParticleSet* ptclA  = &targetPtcl;
  if (sourceInp != targetPtcl.getName())
  {
    //renameProperty(sourceInp);
    auto pit(ptclPool.find(sourceInp));
    if (pit == ptclPool.end())
    {
      ERRORMSG("Missing source ParticleSet" << sourceInp);
      APP_ABORT("HamiltonianFactory::addCoulombPotential");
      return;
    }
    ptclA = (*pit).second;
  }
  if (sourceInp == targetInp) // AA type
  {
    if (!applyPBC && ptclA->getTotalNum() == 1)
    {
      app_log() << "  CoulombAA for " << sourceInp << " is not created.  Number of particles == 1 and nonPeriodic"
                << std::endl;
      return;
    }
    bool quantum = (sourceInp == targetPtcl.getName());
#ifdef QMC_CUDA
    if (applyPBC)
      targetH->addOperator(std::make_unique<CoulombPBCAA_CUDA>(*ptclA, quantum, doForces), title, physical);
    else
    {
      if (quantum)
        targetH->addOperator(std::make_unique<CoulombPotentialAA_CUDA>(*ptclA, true), title, physical);
      else
        targetH->addOperator(std::make_unique<CoulombPotential<Return_t>>(*ptclA, quantum, doForces), title, physical);
    }
#else
    if (applyPBC)
    {
      if (use_gpu.empty())
        use_gpu = ptclA->getCoordinates().getKind() == DynamicCoordinateKind::DC_POS_OFFLOAD ? "yes" : "no";

      if (use_gpu == "yes" && ptclA->getCoordinates().getKind() != DynamicCoordinateKind::DC_POS_OFFLOAD)
        throw std::runtime_error("Requested gpu=yes in CoulombPBCAA but the particle set has gpu=no.");

      targetH->addOperator(std::make_unique<CoulombPBCAA>(*ptclA, quantum, doForces, use_gpu == "yes"), title, physical);
    }
    else
    {
      targetH->addOperator(std::make_unique<CoulombPotential<Return_t>>(*ptclA, quantum, doForces), title, physical);
    }
#endif
  }
  else //X-e type, for X=some other source
  {
#ifdef QMC_CUDA
    if (applyPBC)
      targetH->addOperator(std::make_unique<CoulombPBCAB_CUDA>(*ptclA, targetPtcl), title);
    else
      targetH->addOperator(std::make_unique<CoulombPotentialAB_CUDA>(*ptclA, targetPtcl), title);
#else
    if (applyPBC)
      targetH->addOperator(std::make_unique<CoulombPBCAB>(*ptclA, targetPtcl), title);
    else
      targetH->addOperator(std::make_unique<CoulombPotential<Return_t>>(*ptclA, targetPtcl, true), title);
#endif
  }
}

void HamiltonianFactory::addForceHam(xmlNodePtr cur)
{
#if OHMMS_DIM == 3
  std::string a("ion0"), targetName("e"), title("ForceBase"), pbc("yes"), PsiName = "psi0";
  OhmmsAttributeSet hAttrib;
  std::string mode("bare");
  //hAttrib.add(title,"id");
  hAttrib.add(title, "name");
  hAttrib.add(a, "source");
  hAttrib.add(targetName, "target");
  hAttrib.add(pbc, "pbc");
  hAttrib.add(mode, "mode");
  hAttrib.add(PsiName, "psi");
  hAttrib.put(cur);
  app_log() << "HamFac forceBase mode " << mode << std::endl;
  bool applyPBC = (PBCType && pbc == "yes");

  bool quantum = (a == targetPtcl.getName());

  renameProperty(a);
  auto pit(ptclPool.find(a));
  if (pit == ptclPool.end())
  {
    ERRORMSG("Missing source ParticleSet" << a)
    return;
  }
  ParticleSet* source = (*pit).second;
  pit                 = ptclPool.find(targetName);
  if (pit == ptclPool.end())
  {
    ERRORMSG("Missing target ParticleSet" << targetName)
    return;
  }
  ParticleSet* target = (*pit).second;
  //bool applyPBC= (PBCType && pbc=="yes");
  if (mode == "bare")
  {
    std::unique_ptr<BareForce> bareforce = std::make_unique<BareForce>(*source, *target);
    bareforce->put(cur);
    targetH->addOperator(std::move(bareforce), title, false);
  }
  else if (mode == "cep")
  {
    if (applyPBC == true)
    {
      std::unique_ptr<ForceChiesaPBCAA> force_chi = std::make_unique<ForceChiesaPBCAA>(*source, *target, true);
      force_chi->put(cur);
      targetH->addOperator(std::move(force_chi), title, false);
    }
    else
    {
      std::unique_ptr<ForceCeperley> force_cep = std::make_unique<ForceCeperley>(*source, *target);
      force_cep->put(cur);
      targetH->addOperator(std::move(force_cep), title, false);
    }
  }
  else if (mode == "acforce")
  {
    app_log() << "Adding Assaraf-Caffarel total force.\n";
    auto psi_it(psiPool.find(PsiName));
    if (psi_it == psiPool.end())
    {
      APP_ABORT("Unknown psi \"" + PsiName + "\" for zero-variance force.");
    }
    TrialWaveFunction& psi           = *psi_it->second->getTWF();
    std::unique_ptr<ACForce> acforce = std::make_unique<ACForce>(*source, *target, psi, *targetH);
    acforce->put(cur);
    targetH->addOperator(std::move(acforce), title, false);
  }
  else
  {
    ERRORMSG("Failed to recognize Force mode " << mode);
  }
#endif
}

void HamiltonianFactory::addPseudoPotential(xmlNodePtr cur)
{
#if OHMMS_DIM == 3
  std::string src("i"), title("PseudoPot"), wfname("invalid"), format("xml");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(title, "name");
  pAttrib.add(src, "source");
  pAttrib.add(wfname, "wavefunction");
  pAttrib.add(format, "format"); //temperary tag to switch between format
  pAttrib.put(cur);
  if (format == "old")
  {
    APP_ABORT("pseudopotential Table format is not supported.");
  }
  renameProperty(src);
  renameProperty(wfname);
  auto pit(ptclPool.find(src));
  if (pit == ptclPool.end())
  {
    ERRORMSG("Missing source ParticleSet" << src)
    return;
  }
  ParticleSet* ion = (*pit).second;
  auto oit(psiPool.find(wfname));
  TrialWaveFunction* psi = 0;
  if (oit == psiPool.end())
  {
    if (psiPool.empty())
      return;
    app_error() << "  Cannot find " << wfname << " in the Wavefunction pool. Using the first wavefunction."
                << std::endl;
    psi = (*(psiPool.begin())).second->getTWF();
  }
  else
  {
    psi = (*oit).second->getTWF();
  }
  //remember the TrialWaveFunction used by this pseudopotential
  psiName = wfname;
  app_log() << std::endl << "  ECPotential builder for pseudopotential " << std::endl;
  ECPotentialBuilder ecp(*targetH, *ion, targetPtcl, *psi, myComm);
  ecp.put(cur);
#else
  APP_ABORT("HamiltonianFactory::addPseudoPotential\n pairpot@type=\"pseudo\" is invalid if DIM != 3");
#endif
}

//  void
//  HamiltonianFactory::addConstCoulombPotential(xmlNodePtr cur, std::string& nuclei)
//  {
//    OhmmsAttributeSet hAttrib;
//    std::string hname("IonIon");
//    std::string forces("no");
//    hAttrib.add(forces,"forces");
//    hAttrib.add(hname,"name");
//    hAttrib.put(cur);
//    bool doForces = (forces == "yes") || (forces == "true");
//
//    app_log() << "  Creating Coulomb potential " << nuclei << "-" << nuclei << std::endl;
//    renameProperty(nuclei);
//    PtclPoolType::iterator pit(ptclPool.find(nuclei));
//    if(pit != ptclPool.end()) {
//      ParticleSet* ion=(*pit).second;
//      if(PBCType)
//      {
//#ifdef QMC_CUDA
//	targetH->addOperator(new CoulombPBCAA_CUDA(*ion,false,doForces),hname);
//#else
//	targetH->addOperator(new CoulombPBCAATemp(*ion,false,doForces),hname);
//#endif
//      } else {
//        if(ion->getTotalNum()>1)
//          targetH->addOperator(new CoulombPotential<Return_t>(ion),hname);
//          //targetH->addOperator(new IonIonPotential(*ion),hname);
//      }
//    }
//  }

} // namespace qmcplusplus
