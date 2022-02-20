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
// #include "QMCHamiltonians/ZeroVarObs.h"
#if !defined(QMC_CUDA)
#include "QMCHamiltonians/SkPot.h"
#endif
#include "OhmmsData/AttributeSet.h"
#ifdef QMC_CUDA
#include "QMCHamiltonians/SkEstimator_CUDA.h"
#endif

namespace qmcplusplus
{
HamiltonianFactory::HamiltonianFactory(const std::string& hName,
                                       ParticleSet& qp,
                                       const PtclPoolType& pset,
                                       const PsiPoolType& oset,
                                       Communicate* c)
    : MPIObjectBase(c),
      targetH(std::make_unique<QMCHamiltonian>(hName)),
      targetPtcl(qp),
      ptclPool(pset),
      psiPool(oset),
      myNode(NULL),
      psiName("psi0")
{
  //PBCType is zero or 1 but should be generalized
  PBCType   = targetPtcl.getLattice().SuperCellEnum;
  ClassName = "HamiltonianFactory";
  myName    = hName;
  targetPtcl.set_quantum();
  targetH->addOperator(std::make_unique<BareKineticEnergy>(targetPtcl), "Kinetic");
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
bool HamiltonianFactory::build(xmlNodePtr cur, bool buildtree)
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
  bool attach2Node = false;
  if (buildtree)
  {
    if (myNode == NULL)
    {
      //#if (LIBXMLD_VERSION < 20616)
      //        app_warning() << "   Workaround of libxml2 bug prior to 2.6.x versions" << std::endl;
      //        myNode = xmlCopyNode(cur,2);
      //#else
      //        app_warning() << "   using libxml2 2.6.x versions" << std::endl;
      //        myNode = xmlCopyNode(cur,1);
      //#endif
      myNode = xmlCopyNode(cur, 1);
    }
    else
    {
      attach2Node = true;
    }
  }
  auto psi_it(psiPool.find(psiName));
  if (psi_it == psiPool.end())
    APP_ABORT("Unknown psi \"" + psiName + "\" for target Psi");
  TrialWaveFunction* targetPsi = psi_it->second->getTWF();
  xmlNodePtr cur_saved(cur);
  cur = cur->children;
  while (cur != NULL)
  {
    std::string notype = "0";
    std::string noname = "any";
    std::string cname((const char*)cur->name);
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
    attrib.put(cur);
    renameProperty(sourceInp);
    renameProperty(targetInp);

    int nham = targetH->total_size();
    if (cname == "pairpot")
    {
      if (potType == "coulomb")
      {
        addCoulombPotential(cur);
      }
#if !defined(QMC_CUDA)
      else if (potType == "skpot")
      {
        std::unique_ptr<SkPot> hs = std::make_unique<SkPot>(targetPtcl);
        hs->put(cur);
        targetH->addOperator(std::move(hs), "SkPot", true);
      }
#endif
#if OHMMS_DIM == 3
      else if (potType == "MPC" || potType == "mpc")
        addMPCPotential(cur);
      else if (potType == "pseudo")
        addPseudoPotential(cur);
#endif
    }
    else if (cname == "constant")
    {
      //just to support old input
      if (potType == "coulomb")
        addCoulombPotential(cur);
    }
    else if (cname == "extpot")
    {
      if (potType == "harmonic_ext" || potType == "HarmonicExt")
      {
        std::unique_ptr<HarmonicExternalPotential> hs = std::make_unique<HarmonicExternalPotential>(targetPtcl);
        hs->put(cur);
        targetH->addOperator(std::move(hs), "HarmonicExt", true);
      }
      if (potType == "grid")
      {
        std::unique_ptr<GridExternalPotential> hs = std::make_unique<GridExternalPotential>(targetPtcl);
        hs->put(cur);
        targetH->addOperator(std::move(hs), "Grid", true);
      }
    }
    else if (cname == "estimator")
    {
      if (potType == "flux")
      {
        targetH->addOperator(std::make_unique<ConservedEnergy>(), potName, false);
      }
      else if (potType == "specieskinetic")
      {
        std::unique_ptr<SpeciesKineticEnergy> apot = std::make_unique<SpeciesKineticEnergy>(targetPtcl);
        apot->put(cur);
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
        ParticleSet* target_particle_set = (*pit).second;

        // find source particle set
        auto spit(ptclPool.find(sourceInp));
        if (spit == ptclPool.end())
        {
          APP_ABORT("Unknown source \"" + sourceInp + "\" for LatticeDeviation.");
        }
        ParticleSet* source_particle_set = (*spit).second;

        // read xml node
        OhmmsAttributeSet local_attrib;
        std::string target_group, source_group;
        local_attrib.add(target_group, "tgroup");
        local_attrib.add(source_group, "sgroup");
        local_attrib.put(cur);

        std::unique_ptr<LatticeDeviationEstimator> apot =
            std::make_unique<LatticeDeviationEstimator>(*target_particle_set, *source_particle_set, target_group,
                                                        source_group);
        apot->put(cur);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "Force")
      {
        addForceHam(cur);
      }
      else if (potType == "gofr")
      {
        std::unique_ptr<PairCorrEstimator> apot = std::make_unique<PairCorrEstimator>(targetPtcl, sourceInp);
        apot->put(cur);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "density")
      {
        //          if(PBCType)//only if perioidic
        {
          std::unique_ptr<DensityEstimator> apot = std::make_unique<DensityEstimator>(targetPtcl);
          apot->put(cur);
          targetH->addOperator(std::move(apot), potName, false);
        }
      }
      else if (potType == "spindensity")
      {
        app_log() << "  Adding SpinDensity" << std::endl;
        std::unique_ptr<SpinDensity> apot = std::make_unique<SpinDensity>(targetPtcl);
        apot->put(cur);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "structurefactor")
      {
        app_log() << "  Adding StaticStructureFactor" << std::endl;
        std::unique_ptr<StaticStructureFactor> apot = std::make_unique<StaticStructureFactor>(targetPtcl);
        apot->put(cur);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "orbitalimages")
      {
        app_log() << "  Adding OrbitalImages" << std::endl;
        std::unique_ptr<OrbitalImages> apot =
            std::make_unique<OrbitalImages>(targetPtcl, ptclPool, myComm, *psi_it->second);
        apot->put(cur);
        targetH->addOperator(std::move(apot), potName, false);
      }
#if !defined(REMOVE_TRACEMANAGER)
      else if (potType == "energydensity" || potType == "EnergyDensity")
      {
        app_log() << "  Adding EnergyDensityEstimator" << std::endl;
        std::unique_ptr<EnergyDensityEstimator> apot = std::make_unique<EnergyDensityEstimator>(ptclPool, defaultKE);
        apot->put(cur);
        targetH->addOperator(std::move(apot), potName, false);
      }
      else if (potType == "dm1b")
      {
        app_log() << "  Adding DensityMatrices1B" << std::endl;
        std::string source = "";
        OhmmsAttributeSet attrib;
        attrib.add(source, "source");
        attrib.put(cur);
        auto pit(ptclPool.find(source));
        ParticleSet* Pc = nullptr;
        if (source == "")
          Pc = nullptr;
        else if (pit != ptclPool.end())
          Pc = pit->second;
        else
        {
          APP_ABORT("Unknown source \"" + source + "\" for DensityMatrices1B");
        }
        std::unique_ptr<DensityMatrices1B> apot =
            std::make_unique<DensityMatrices1B>(targetPtcl, *targetPsi, Pc, *psi_it->second);
        apot->put(cur);
        targetH->addOperator(std::move(apot), potName, false);
      }
#endif
      else if (potType == "sk")
      {
        if (PBCType) //only if perioidic
        {
#ifdef QMC_CUDA
          std::unique_ptr<SkEstimator_CUDA> apot = std::make_unique<SkEstimator_CUDA>(targetPtcl);
#else
          std::unique_ptr<SkEstimator> apot = std::make_unique<SkEstimator>(targetPtcl);
#endif
          apot->put(cur);
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
        hAttrib.put(cur);
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
        const TrialWaveFunction& psi             = *psi_it->second->getTWF();
        std::unique_ptr<ChiesaCorrection> chiesa = std::make_unique<ChiesaCorrection>(source, psi);
        targetH->addOperator(std::move(chiesa), "KEcorr", false);
      }
      else if (potType == "skall")
      {
        std::string SourceName = "";
        OhmmsAttributeSet attrib;
        attrib.add(SourceName, "source");
        attrib.put(cur);

        auto pit(ptclPool.find(SourceName));
        if (pit == ptclPool.end())
        {
          APP_ABORT("Unknown source \"" + SourceName + "\" for SkAll.");
        }
        ParticleSet* source = (*pit).second;

        if (PBCType)
        {
          std::unique_ptr<SkAllEstimator> apot = std::make_unique<SkAllEstimator>(*source, targetPtcl);
          apot->put(cur);
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
          BP->put(cur);
          targetH->addOperator(std::move(BP), "Pressure", false);
          int nlen(100);
          attrib.add(nlen, "truncateSum");
          attrib.put(cur);
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
        hAttrib.put(cur);
        auto psi_it(psiPool.find(PsiName));
        if (psi_it == psiPool.end())
        {
          APP_ABORT("Unknown psi \"" + PsiName + "\" for momentum.");
        }
        TrialWaveFunction* psi                = (*psi_it).second->getTWF();
        std::unique_ptr<MomentumEstimator> ME = std::make_unique<MomentumEstimator>(targetPtcl, *psi);
        bool rt(myComm->rank() == 0);
        ME->putSpecial(cur, targetPtcl, rt);
        targetH->addOperator(std::move(ME), "MomentumEstimator", false);
      }
    }
    else if (cname == "Kinetic")
    {
      std::string TargetName = "e";
      std::string SourceName = "I";
      OhmmsAttributeSet hAttrib;
      hAttrib.add(TargetName, "Dependant");
      hAttrib.add(SourceName, "Independent");
      hAttrib.put(cur);
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

    if (attach2Node)
      xmlAddChild(myNode, xmlCopyNode(cur, 1));
    cur = cur->next;
  }
  //add observables with physical and simple estimators
  targetH->addObservables(targetPtcl);
  //do correction
  bool dmc_correction = false;
  cur                 = cur_saved->children;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    std::string potType("0");
    OhmmsAttributeSet attrib;
    attrib.add(potType, "type");
    attrib.put(cur);
    if (cname == "estimator")
    {
      if (potType == "ZeroVarObs")
      {
        app_log() << "  Not Adding ZeroVarObs Operator" << std::endl;
        //         ZeroVarObs* FW=new ZeroVarObs();
        //         FW->put(cur,*targetH,targetPtcl);
        //         targetH->addOperator(FW,"ZeroVarObs",false);
      }
      //         else if(potType == "DMCCorrection")
      //         {
      //           TrialDMCCorrection* TE = new TrialDMCCorrection();
      //           TE->putSpecial(cur,*targetH,targetPtcl);
      //           targetH->addOperator(TE,"DMC_CORR",false);
      //           dmc_correction=true;
      //         }
      else if (potType == "ForwardWalking")
      {
        app_log() << "  Adding Forward Walking Operator" << std::endl;
        std::unique_ptr<ForwardWalking> FW = std::make_unique<ForwardWalking>();
        FW->putSpecial(cur, *targetH, targetPtcl);
        targetH->addOperator(std::move(FW), "ForwardWalking", false);
        dmc_correction = true;
      }
    }
    cur = cur->next;
  }
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
  bool success = build(cur, false);
  return success;
}

} // namespace qmcplusplus
