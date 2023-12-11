//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Bryan Clark, bclark@Princeton.edu, Princeton University
//                    D.C. Yang, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file HamiltonianFactory.cpp
 *@brief Definition of a HamiltonianFactory
 */
#include "HamiltonianFactory.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "QMCHamiltonians/ConservedEnergy.h"
#include "QMCHamiltonians/SpeciesKineticEnergy.h"
#include "QMCHamiltonians/LatticeDeviationEstimator.h"
#include "QMCHamiltonians/MomentumEstimator.h"
#include "QMCHamiltonians/Pressure.h"
#include "QMCHamiltonians/ForwardWalking.h"
#include "QMCHamiltonians/PairCorrEstimator.h"
#include "QMCHamiltonians/DensityEstimator.h"
#include "QMCHamiltonians/SkEstimator.h"
#include "QMCHamiltonians/HarmonicExternalPotential.h"
#include "QMCHamiltonians/GridExternalPotential.h"
#include "QMCHamiltonians/StaticStructureFactor.h"
#include "QMCHamiltonians/SpinDensity.h"
#include "QMCHamiltonians/OrbitalImages.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "QMCHamiltonians/EnergyDensityEstimator.h"
#include "QMCHamiltonians/DensityMatrices1B.h"
#endif
#if OHMMS_DIM == 3
#include "QMCHamiltonians/ChiesaCorrection.h"
#include "QMCHamiltonians/SkAllEstimator.h"
#endif
#include "QMCHamiltonians/SkPot.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
HamiltonianFactory::HamiltonianFactory(const std::string& hName,
                                       ParticleSet& qp,
                                       const PSetMap& pset,
                                       const PsiPoolType& oset,
                                       Communicate* c)
    : MPIObjectBase(c),
      targetH(std::make_unique<QMCHamiltonian>(hName)),
      targetPtcl(qp),
      ptclPool(pset),
      psiPool(oset),
      psiName("psi0")
{
  //PBCType is zero or 1 but should be generalized
  PBCType   = targetPtcl.getLattice().SuperCellEnum;
  ClassName = "HamiltonianFactory";
  myName    = hName;
  targetPtcl.set_quantum();
}

/** main hamiltonian build function
 * @param cur element node <hamiltonian/>
 * @param buildtree if true, build xml tree for a reuse
 *
 * A valid hamiltonian node contains
 * \xmlonly
 *  <hamiltonian target="e">
 *    <pairpot type="coulomb" name="ElecElec" source="e"/>
 *    <pairpot type="coulomb" name="IonElec" source="i"/>
 *    <pairpot type="coulomb" name="IonIon" source="i" target="i"/>
 *  </hamiltonian>
 * \endxmlonly
 */
bool HamiltonianFactory::build(xmlNodePtr cur)
{
  if (cur == NULL)
    return false;

  app_summary() << std::endl;
  app_summary() << " Hamiltonian and observables" << std::endl;
  app_summary() << " ---------------------------" << std::endl;
  app_summary() << "  Name: " << myName << std::endl;
  app_summary() << std::endl;

  std::string htype("generic"), source("i"), defaultKE("yes");
  OhmmsAttributeSet hAttrib;
  hAttrib.add(htype, "type");
  hAttrib.add(source, "source");
  hAttrib.add(defaultKE, "default");
  hAttrib.put(cur);
  renameProperty(source);
  auto psi_it(psiPool.find(psiName));
  if (psi_it == psiPool.end())
    APP_ABORT("Unknown psi \"" + psiName + "\" for target Psi");
  TrialWaveFunction* targetPsi = psi_it->second.get();
  // KineticEnergy must be the first element in the hamiltonian array.
  if (defaultKE != "no")
    targetH->addOperator(std::make_unique<BareKineticEnergy>(targetPtcl, *targetPsi), "Kinetic");

  // Virtual particle sets only need to carry distance tables used by the wavefunction.
  // Other Hamiltonian elements or estimators may add distance tables in particle sets.
  // Process pseudopotentials first to minimize the needed distance tables in virtual particle sets.
  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "pairpot" && getXMLAttributeValue(element, "type") == "pseudo")
      addPseudoPotential(element);
  });

  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    std::string notype = "0";
    std::string noname = "any";
    std::string potType(notype);
    std::string potName(noname);
    std::string potUnit("hartree");
    std::string estType("coulomb");
    std::string sourceInp(targetPtcl.getName());
    std::string targetInp(targetPtcl.getName());
    OhmmsAttributeSet attrib;
    attrib.add(sourceInp, "source");
    attrib.add(sourceInp, "sources");
    attrib.add(targetInp, "target");
    attrib.add(potType, "type");
    attrib.add(potName, "name");
    attrib.add(potUnit, "units");
    attrib.add(estType, "potential");
    attrib.put(element);
    renameProperty(sourceInp);
    renameProperty(targetInp);

    int nham = targetH->total_size();
    if (cname == "pairpot")
    {
      if (potType == "coulomb")
        addCoulombPotential(element);
      else if (potType == "skpot")
      {
        std::unique_ptr<SkPot> hs = std::make_unique<SkPot>(targetPtcl);
        hs->put(element);
        targetH->addOperator(std::move(hs), "SkPot", true);
      }
#if OHMMS_DIM == 3
      else if (potType == "MPC" || potType == "mpc")
        addMPCPotential(element);
#endif
    }
    else if (cname == "constant")
    {
      //just to support old input
      if (potType == "coulomb")
        addCoulombPotential(element);
    }
    else if (cname == "extpot")
    {
      if (potType == "harmonic_ext" || potType == "HarmonicExt")
      {
        std::unique_ptr<HarmonicExternalPotential> hs = std::make_unique<HarmonicExternalPotential>(targetPtcl);
        hs->put(element);
        targetH->addOperator(std::move(hs), "HarmonicExt", true);
      }
      if (potType == "grid")
      {
        std::unique_ptr<GridExternalPotential> hs = std::make_unique<GridExternalPotential>(targetPtcl);
        hs->put(element);
        targetH->addOperator(std::move(hs), "Grid", true);
      }
    }
    else if (cname == "estimator")
    {
      if (potType == "flux")
        targetH->addOperator(std::make_unique<ConservedEnergy>(), potName, false);
      else if (potType == "specieskinetic")
      {
        std::unique_ptr<SpeciesKineticEnergy> apot = std::make_unique<SpeciesKineticEnergy>(targetPtcl);
        apot->put(element);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "latticedeviation")
      {
        // find target particle set
        auto pit(ptclPool.find(targetInp));
        if (pit == ptclPool.end())
        {
          APP_ABORT("Unknown target \"" + targetInp + "\" for LatticeDeviation.");
        }

        // find source particle set
        auto spit(ptclPool.find(sourceInp));
        if (spit == ptclPool.end())
        {
          APP_ABORT("Unknown source \"" + sourceInp + "\" for LatticeDeviation.");
        }

        // read xml node
        OhmmsAttributeSet local_attrib;
        std::string target_group, source_group;
        local_attrib.add(target_group, "tgroup");
        local_attrib.add(source_group, "sgroup");
        local_attrib.put(element);

        std::unique_ptr<LatticeDeviationEstimator> apot =
            std::make_unique<LatticeDeviationEstimator>(*pit->second, *spit->second, target_group, source_group);
        apot->put(element);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "Force")
        addForceHam(element);
      else if (potType == "gofr")
      {
        std::unique_ptr<PairCorrEstimator> apot = std::make_unique<PairCorrEstimator>(targetPtcl, sourceInp);
        apot->put(element);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "density")
      {
        std::unique_ptr<DensityEstimator> apot = std::make_unique<DensityEstimator>(targetPtcl);
        apot->put(element);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "spindensity")
      {
        app_log() << "  Adding SpinDensity" << std::endl;
        std::unique_ptr<SpinDensity> apot = std::make_unique<SpinDensity>(targetPtcl);
        apot->put(element);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "structurefactor")
      {
        app_log() << "  Adding StaticStructureFactor" << std::endl;
        std::unique_ptr<StaticStructureFactor> apot = std::make_unique<StaticStructureFactor>(targetPtcl);
        apot->put(element);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "orbitalimages")
      {
        app_log() << "  Adding OrbitalImages" << std::endl;
        std::unique_ptr<OrbitalImages> apot =
            std::make_unique<OrbitalImages>(targetPtcl, ptclPool, myComm, targetPsi->getSPOMap());
        apot->put(element);
        targetH->addOperator(std::move(apot), potName, false);
      }
#if !defined(REMOVE_TRACEMANAGER)
      else if (potType == "energydensity" || potType == "EnergyDensity")
      {
        app_log() << "  Adding EnergyDensityEstimator" << std::endl;
        std::unique_ptr<EnergyDensityEstimator> apot = std::make_unique<EnergyDensityEstimator>(ptclPool, defaultKE);
        apot->put(element);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "dm1b")
      {
        app_log() << "  Adding DensityMatrices1B" << std::endl;
        std::string source = "";
        OhmmsAttributeSet attrib;
        attrib.add(source, "source");
        attrib.put(element);
        auto pit(ptclPool.find(source));
        ParticleSet* Pc = nullptr;
        if (source == "")
          Pc = nullptr;
        else if (pit != ptclPool.end())
          Pc = pit->second.get();
        else
        {
          APP_ABORT("Unknown source \"" + source + "\" for DensityMatrices1B");
        }
        std::unique_ptr<DensityMatrices1B> apot = std::make_unique<DensityMatrices1B>(targetPtcl, *targetPsi, Pc);
        apot->put(element);
        targetH->addOperator(std::move(apot), potName, false);
      }
#endif
      else if (potType == "sk")
      {
        if (PBCType) //only if perioidic
        {
          std::unique_ptr<SkEstimator> apot = std::make_unique<SkEstimator>(targetPtcl);
          apot->put(element);
          targetH->addOperator(std::move(apot), potName, false);
          app_log() << "Adding S(k) estimator" << std::endl;
        }
      }
#if OHMMS_DIM == 3
      else if (potType == "chiesa")
      {
        std::string PsiName    = "psi0";
        std::string SourceName = "e";
        OhmmsAttributeSet hAttrib;
        hAttrib.add(PsiName, "psi");
        hAttrib.add(SourceName, "source");
        hAttrib.put(element);
        auto pit(ptclPool.find(SourceName));
        if (pit == ptclPool.end())
        {
          APP_ABORT("Unknown source \"" + SourceName + "\" for Chiesa correction.");
        }
        ParticleSet& source = *pit->second;
        auto psi_it(psiPool.find(PsiName));
        if (psi_it == psiPool.end())
        {
          APP_ABORT("Unknown psi \"" + PsiName + "\" for Chiesa correction.");
        }
        const TrialWaveFunction& psi             = *psi_it->second;
        std::unique_ptr<ChiesaCorrection> chiesa = std::make_unique<ChiesaCorrection>(source, psi);
        targetH->addOperator(std::move(chiesa), "KEcorr", false);
      }
      else if (potType == "skall")
      {
        std::string SourceName = "";
        OhmmsAttributeSet attrib;
        attrib.add(SourceName, "source");
        attrib.put(element);

        auto pit(ptclPool.find(SourceName));
        if (pit == ptclPool.end())
        {
          APP_ABORT("Unknown source \"" + SourceName + "\" for SkAll.");
        }

        if (PBCType)
        {
          std::unique_ptr<SkAllEstimator> apot = std::make_unique<SkAllEstimator>(*pit->second, targetPtcl);
          apot->put(element);
          targetH->addOperator(std::move(apot), potName, false);
          app_log() << "Adding S(k) ALL estimator" << std::endl;
        }
      }

#endif
      else if (potType == "Pressure")
      {
        if (estType == "coulomb")
        {
          std::unique_ptr<Pressure> BP = std::make_unique<Pressure>(targetPtcl);
          BP->put(element);
          targetH->addOperator(std::move(BP), "Pressure", false);
          int nlen(100);
          attrib.add(nlen, "truncateSum");
          attrib.put(element);
          //             DMCPressureCorr* DMCP = new DMCPressureCorr(targetPtcl,nlen);
          //             targetH->addOperator(DMCP,"PressureSum",false);
        }
      }
      else if (potType == "momentum")
      {
        app_log() << "  Adding Momentum Estimator" << std::endl;
        std::string PsiName = "psi0";
        OhmmsAttributeSet hAttrib;
        hAttrib.add(PsiName, "wavefunction");
        hAttrib.put(element);
        auto psi_it(psiPool.find(PsiName));
        if (psi_it == psiPool.end())
        {
          APP_ABORT("Unknown psi \"" + PsiName + "\" for momentum.");
        }
        std::unique_ptr<MomentumEstimator> ME = std::make_unique<MomentumEstimator>(targetPtcl, *psi_it->second);
        bool rt(myComm->rank() == 0);
        ME->putSpecial(element, targetPtcl, rt);
        targetH->addOperator(std::move(ME), "MomentumEstimator", false);
      }
    }

    if (nham < targetH->total_size()) //if(cname!="text" && cname !="comment")
    {
      if (potName == noname)
      {
        potName = potType;
        app_log() << "Provide name for hamiltonian element for type " << potType << std::endl;
      }
      //APP_ABORT("HamiltonianFactory::build\n  a name for operator of type "+cname+" "+potType+" must be provided in the xml input");
      targetH->addOperatorType(potName, potType);
    }
  });

  //add observables with physical and simple estimators
  targetH->addObservables(targetPtcl);
  //do correction
  bool dmc_correction = false;
  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    std::string potType("0");
    OhmmsAttributeSet attrib;
    attrib.add(potType, "type");
    attrib.put(element);
    if (cname == "estimator" && potType == "ForwardWalking")
    {
      app_log() << "  Adding Forward Walking Operator" << std::endl;
      std::unique_ptr<ForwardWalking> FW = std::make_unique<ForwardWalking>();
      FW->putSpecial(element, *targetH, targetPtcl);
      targetH->addOperator(std::move(FW), "ForwardWalking", false);
      dmc_correction = true;
    }
  });
  //evaluate the observables again
  if (dmc_correction)
    targetH->addObservables(targetPtcl);
  return true;
}


void HamiltonianFactory::renameProperty(const std::string& a, const std::string& b) { RenamedProperty[a] = b; }

void HamiltonianFactory::renameProperty(std::string& aname)
{
  std::map<std::string, std::string>::iterator it(RenamedProperty.find(aname));
  if (it != RenamedProperty.end())
  {
    aname = (*it).second;
  }
}

bool HamiltonianFactory::put(xmlNodePtr cur)
{
  bool success = build(cur);
  return success;
}

} // namespace qmcplusplus
