//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file WaveFunctionFactory.cpp
 *@brief Definition of a WaveFunctionFactory
 */
#include "WaveFunctionFactory.h"
#include "QMCWaveFunctions/Jastrow/JastrowBuilder.h"
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "QMCWaveFunctions/LatticeGaussianProductBuilder.h"
#include "QMCWaveFunctions/ExampleHeBuilder.h"

#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#else
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#endif

#include "QMCWaveFunctions/PlaneWave/PWOrbitalBuilder.h"
#if OHMMS_DIM == 3 && !defined(QMC_COMPLEX)
#include "QMCWaveFunctions/AGPDeterminantBuilder.h"
#endif

#include "Utilities/ProgressReportEngine.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"
namespace qmcplusplus
{
WaveFunctionFactory::WaveFunctionFactory(ParticleSet& qp,
                                         const PSetMap& pset,
                                         Communicate* c)
    : MPIObjectBase(c),
      targetPtcl(qp),
      ptclPool(pset)
{
  ClassName = "WaveFunctionFactory";
}

WaveFunctionFactory::~WaveFunctionFactory() = default;

std::unique_ptr<TrialWaveFunction> WaveFunctionFactory::buildTWF(xmlNodePtr cur)
{
  // YL: how can this happen?
  if (cur == NULL)
    return nullptr;

  ReportEngine PRE(ClassName, "build");

  std::string psiName("psi0"), tasking;
  OhmmsAttributeSet pAttrib;
  pAttrib.add(psiName, "id");
  pAttrib.add(psiName, "name");
  pAttrib.add(tasking, "tasking", {"no", "yes"});
  pAttrib.put(cur);

  app_summary() << std::endl;
  app_summary() << " Many-body wavefunction" << std::endl;
  app_summary() << " -------------------" << std::endl;
  app_summary() << "  Name: " << psiName << "   Tasking: " << (tasking == "yes" ? "yes" : "no") << std::endl;
  app_summary() << std::endl;

  auto targetPsi = std::make_unique<TrialWaveFunction>(psiName, tasking == "yes");
  targetPsi->setMassTerm(targetPtcl);
  targetPsi->storeXMLNode(cur);

  SPOSetBuilderFactory sposet_builder_factory(myComm, targetPtcl, ptclPool);

  std::string vp_file_to_load;
  cur          = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "sposet_builder" || cname == "sposet_collection")
      sposet_builder_factory.buildSPOSetCollection(cur);
    else if (cname == WaveFunctionComponentBuilder::detset_tag)
    {
      addFermionTerm(*targetPsi, sposet_builder_factory, cur);
      bool foundtwist(false);
      xmlNodePtr kcur = cur->children;
      while (kcur != NULL)
      {
        std::string kname((const char*)(kcur->name));
        if (kname == "h5tag")
        {
          std::string hdfName;
          OhmmsAttributeSet attribs;
          attribs.add(hdfName, "name");
          if (hdfName == "twistAngle")
          {
            std::vector<ParticleSet::RealType> tsts(3, 0);
            putContent(tsts, kcur);
            targetPsi->setTwist(tsts);
            foundtwist = true;
          }
        }
        kcur = kcur->next;
      }
      if (!foundtwist)
      {
        //default twist is [0 0 0]
        std::vector<ParticleSet::RealType> tsts(3, 0);
        targetPsi->setTwist(tsts);
      }
    }
    else if (cname == WaveFunctionComponentBuilder::jastrow_tag)
    {
      auto jbuilder = std::make_unique<JastrowBuilder>(myComm, targetPtcl, ptclPool);
      targetPsi->addComponent(jbuilder->buildComponent(cur));
    }
    else if (cname == "fdlrwfn")
    {
      APP_ABORT("FDLR wave functions are not currently supported.");
    }
    else if (cname == WaveFunctionComponentBuilder::ionorb_tag)
    {
      auto builder = std::make_unique<LatticeGaussianProductBuilder>(myComm, targetPtcl, ptclPool);
      targetPsi->addComponent(builder->buildComponent(cur));
    }
    else if ((cname == "Molecular") || (cname == "molecular"))
    {
      APP_ABORT("  Removed Helium Molecular terms from qmcpack ");
    }
    else if (cname == "example_he")
    {
      auto exampleHe_builder = std::make_unique<ExampleHeBuilder>(myComm, targetPtcl, ptclPool);
      targetPsi->addComponent(exampleHe_builder->buildComponent(cur));
    }
#if !defined(QMC_COMPLEX) && OHMMS_DIM == 3
    else if (cname == "agp")
    {
      auto agpbuilder = std::make_unique<AGPDeterminantBuilder>(myComm, targetPtcl, ptclPool);
      targetPsi->addComponent(agpbuilder->buildComponent(cur));
    }
#endif
    else if (cname == "override_variational_parameters")
    {
      OhmmsAttributeSet attribs;
      attribs.add(vp_file_to_load, "href");
      attribs.put(cur);
    }

    cur = cur->next;
  }
  //{
  //  ReportEngine PREA("TrialWaveFunction","print");
  //  targetPsi->VarList.print(app_log());
  //}
  // synch all parameters. You don't want some being different if same name.
  opt_variables_type dummy;
  targetPsi->checkInVariables(dummy);
  dummy.resetIndex();
  targetPsi->checkOutVariables(dummy);

  if (!vp_file_to_load.empty())
  {
    app_log() << "  Reading variational parameters from " << vp_file_to_load << std::endl;
    dummy.readFromHDF(vp_file_to_load);
  }

  targetPsi->resetParameters(dummy);
  targetPsi->storeSPOMap(sposet_builder_factory.exportSPOSets());
  return targetPsi;
}

bool WaveFunctionFactory::addFermionTerm(TrialWaveFunction& psi, SPOSetBuilderFactory& spo_factory, xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "addFermionTerm");
  std::string orbtype("MolecularOrbital");
  std::string nuclei("i");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(orbtype, "type");
  oAttrib.add(nuclei, "source");
  oAttrib.put(cur);
  std::unique_ptr<WaveFunctionComponentBuilder> detbuilder;
  if (orbtype == "electron-gas")
  {
#if defined(QMC_COMPLEX)
    detbuilder = std::make_unique<ElectronGasComplexOrbitalBuilder>(myComm, targetPtcl);
#else
    detbuilder = std::make_unique<ElectronGasOrbitalBuilder>(myComm, targetPtcl);
#endif
  }
  else if (orbtype == "PWBasis" || orbtype == "PW" || orbtype == "pw")
  {
    detbuilder = std::make_unique<PWOrbitalBuilder>(myComm, targetPtcl, ptclPool);
  }
  else
    detbuilder = std::make_unique<SlaterDetBuilder>(myComm, spo_factory, targetPtcl, psi, ptclPool);
  psi.addComponent(detbuilder->buildComponent(cur));
  return true;
}

} // namespace qmcplusplus
